import argparse
import re
import sys,os
import subprocess
import traceback

from collections.abc import Mapping
from copy import copy, deepcopy
from datetime import datetime
from fractions import Fraction
from itertools import chain
from pathlib import Path
from textwrap import TextWrapper, indent, dedent, fill

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

from permute import Permutation, Span, Sn, Trivial, Distributions


import logging
logger = logging.getLogger("ChPT")
logextra = {'where': ''}
VERBOSE = logging.INFO - 1
logging.addLevelName(VERBOSE, "VERBOSE")
# From https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
# ... with some customization
class ColorFormatter(logging.Formatter):

    def __init__(self, format):
        super().__init__()

        BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
        RESET = "\033[0m"
        COLOR = "\033[3%dm"
        BOLD = "\033[1m"
        REVERSE = "\033[7m"

        self.formats = {
            logging.DEBUG: COLOR%BLUE + format + RESET,
                    VERBOSE: COLOR%WHITE + format + RESET,
            logging.INFO: COLOR%WHITE + BOLD + format + RESET,
            logging.WARNING: COLOR%YELLOW + BOLD + format + RESET,
            logging.ERROR: COLOR%RED + BOLD + format + RESET,
            logging.CRITICAL: COLOR%RED + REVERSE + format + RESET,
            None: format
        }

    def format(self, record):
        return logging.Formatter(fmt=self.formats.get(record.levelno, self.formats[None])).format(record)

AUTHOR = "ChPTlib, by Mattias Sjö 2025"

PARTICLES = {
    'M': 'meson',
    'V': 'vector',
    'A': 'axial vector',
    'S': 'scalar',
    'P': 'pseudoscalar'}

# Since \n is not allowed in f-strings, use this instead.
# Supports various indentation options as well for convenience.
def newline(indent=0, *, prefix='', infix='', suffix='', indent_str=' '*4):
    return prefix + '\n' + infix + indent_str*indent + suffix
# Like ', '.join(items) but with final "and" (Oxford comma by default)
def and_join(items, oxford=True):
    items = list(items)
    match len(items):
        case 0:
            return ""
        case 1:
            return items[0]
        case 2:
            return f"{items[0]} and {items[1]}"
        case _:
            return f"{', '.join(items[:-1])}{',' if oxford else ''} and {items[-1]}"

class ChPTError(Exception):
    pass

# Convert power-counting orders between O(p^X) format and NXLO format
def order_OtoN(order):
    if order < 2 or order % 2:
        raise ChPTError(f"Invalid power-counting order: O(p^{order})")
    if order > 6:
        return f"N{order//2-1}LO"
    else:
        return f"{'N'*(order//2-1)}LO"
def order_NtoO(order):
    if re.fullmatch('N*LO', order):
        return 2*len(order) - 2
    elif re.fullmatch('N[0-9]+LO', order):
        return 2*int(order[1:-2]) + 2
    else:
        raise ChPTError(f"Invalid power-counting order: '{order}'")

class Diagram:
    def __init__(self, name, n_ext):

        logger.debug(f"Initializing {n_ext}-leg diagram {name}", extra=logextra)
        self.name = name
        self.n_legs = n_ext
        self.symmetry_factor = {None: '1'}
        self.permutation_group = None
        self.flags = []

        self.vertices = []
        self.vertex_name_map = {}
        self.edges = []
        self.external = []

    def get_vertex(self, name_or_index):
        if string in self.vertex_name_map:
            return self.vertices[ self.vertex_name_map[name_or_index]  ]
        else:
            return self.vertices[int(name_or_index) - 1]

    def add_vertex(self, vertex):
        if vertex.name:
            if vertex.name in self.vertex_name_map:
                raise ChPTError(f"Ambiguous vertex name: '{vertex.name}'")
            self.vertex_name_map[vertex.name] = len(self.vertices)

        self.vertices.append(vertex)
        logger.debug(f"Added vertex {len(self.vertices)}{f' ({vertex.name})' if vertex.name else ''}: {vertex}", extra=logextra)

        # This should not carry over when copying a vertex,
        # and it's easier to zero it here than to mess with copy()
        vertex.ingoing_momenta = []

    def add_edge(self, edge):
        # Edge constructor performs bounds checks but imports can circumvent them
        if edge.source >= len(self.vertices) or edge.destination >= len(self.vertices):
            raise ChPTError(f"Importing edge {edge} into {len(self.vertices)}-vertex diagram")
        self.edges.append(edge)
        logger.debug(f"Added edge {edge}", extra=logextra)

    def set_symmetry_factor(self, tokens):
        self.symmetry_factor[None] = tokens[0];
        for token in tokens[1:]:
            split = token.split(':')
            if len(split) != 2 or not split[0] or not split[1]:
                raise ChPTError("Alternative symmetry factor should be of the format VARIABLE:FACTOR")
            self.symmetry_factor[split[0]] = split[1]

    def set_permutation_group(self, tokens):
        if self.permutation_group is not None:
            raise ChPTError("Multiple symmetry group specifications in the same diagram")
        if not tokens or not tokens[0]:
            raise ChPTError("Empty symmetry group specification")

        if tokens[0][0] == '[':
            if len(tokens) != 1:
                raise ChPTError("Only a single distribution may be specified")
            if tokens[0][-1] != ']' or '[' in tokens[0][1:] or ']' in tokens[0][:-1]:
                raise ChPTError(f"Sincle bracket-enclosed list expected, got '{tokens[0]}'")

            sets = []
            equiv = []
            for token in tokens[0][1:-1].split(','):
                if not token or not token[0].isdigit():
                    raise ChPTError(f"Number (and possible label) expected, got '{tokens}'")

                index = 1
                while index < len(token) and token[index].isdigit():
                    index += 1
                sets.append(int(token[:index]))

                if index >= len(token):
                    equiv.append(True)
                elif token[index:] == '*':
                    equiv.append(False)
                elif token[index:].isalnum():
                    equiv.append(token[index:])
                else:
                    raise ChPTError(f"Malformed set label: expected '*' or alphanumeric string, got '{token[index:]}'")

            if sum(sets) != self.n_legs:
                raise ChPTError(f"Distribution should distribute {self.n_legs} external momenta, got {sum(sets)}")

            self.permutation_group = Distributions(sets, equiv)

        else:
            self.permutation_group = Span(
                Permutation.parse_cycles(size=self.n_legs, string=token, sep=',', base=1) for token in tokens
                ).permutations()

    def add_flags(self, tokens):
        for token in tokens:
            eq = token.find('=')
            if eq < 0:
                self.flags.append((token, None))
            else:
                self.flags.append((token[:eq], token[eq+1:]))

    def legs_up_to(self, legs, particle):
        count = 0
        for p in PARTICLES:
            count += legs.get(p, 0)
            if p == particle:
                return count

    def connect_edge(self, edge, placeholder=False):
        if placeholder:
            self.vertices[edge.source]     .ingoing_momenta.append(None)
            self.vertices[edge.destination].ingoing_momenta.append(None)
        else:
            self.vertices[edge.source]     .ingoing_momenta.append(
                Momentum(f"-({edge.propagator.momentum})", implicit=True))
            self.vertices[edge.destination].ingoing_momenta.append(
                Momentum(f"+({edge.propagator.momentum})", implicit=True))

    def connect_momenta(self, propagators, ext_legs, momenta):
        # Clear things in case this diagram is partly copied from another
        self.external.clear()
        for vertex in self.vertices:
            vertex.ingoing_momenta.clear()

        needs_resolution = False

        for edge in self.edges:
            if edge.unresolved:
                needs_resolution = True
                self.connect_edge(edge, placeholder=True)
            else:
                self.connect_edge(edge)

        for particle in PARTICLES:
            for i, vertex in enumerate(self.vertices):
                count = self.legs_up_to(vertex.legs, particle)
                if count < len(vertex.ingoing_momenta):
                    raise ChPTError(f"Mismatched number of external {PARTICLES[particle]} legs in diagram '{self.name}' (discovered at vertex '{i+1}')")
                for _ in range(count - len(vertex.ingoing_momenta)):
                    self.external.append(i);
                    vertex.ingoing_momenta.append(Momentum(f"+p{len(self.external)}", implicit=True))
                    logger.debug(f"Connected external leg {len(self.external)} ({PARTICLES[particle]}) to vertex {i+1}", extra=logextra)
            if len(self.external) != self.legs_up_to(ext_legs, particle):
                raise ChPTError(f"Mismatched number of external {PARTICLES[particle]} legs in diagram '{self.name}'")

        if needs_resolution:
            self.resolve_momenta(momenta)

    def resolve_momenta(self, momenta):

        # Remove placeholders for unresolved propagators
        # These were used to get the number of external momenta right
        # but are just in the way now
        for vertex in self.vertices:
            vertex.ingoing_momenta = [q for q in vertex.ingoing_momenta if q is not None]

        while True:
            more_needed = False
            made_progress = False

            for edge in self.edges:
                if not edge.unresolved:
                    continue

                # Looking only at source vertices for simplicity
                # There should be no scenario where this hinders resolution
                vertex = self.vertices[edge.source]
                deficit =  vertex.n_legs() - len(vertex.ingoing_momenta)
                if deficit < 0:
                    raise ChPTError(f"{-deficit} too many momenta in vertex {vertex}")
                elif deficit > 1:
                    more_needed = True
                elif deficit == 1:
                    made_progress = True

                    edge.propagator = Propagator(
                            [sum(vertex.ingoing_momenta), edge.propagator],
                            momenta)
                    edge.unresolved = False
                    logger.debug(f"Resolved momentum routing for edge {edge}", extra=logextra)

                    self.connect_edge(edge)

            if not more_needed:
                break
            if not made_progress:
                raise ChPTError(f"Failed to resolve momentum routing in {self.name}")


    def check_CoM(self, propagators, symbols, replacements, shortcircuit=False):
        locals().update(symbols)
        violations = []
        for perm in self.permutation_group:
            for i, vertex in enumerate(self.vertices):
                zero = sum(vertex.ingoing_momenta)
                logging.debug(f"CoM, diagram {self.name}, vertex {i+1}, permutation {perm}: {zero}", extra=logextra)
                zero.substitute({f'p{j+1}' : f'p{perm[j]+1}' for j in range(len(perm))})
                zero.substitute(replacements)
                if zero != 0:
                    if shortcircuit:
                        return False
                    violations.append( (i+1, perm, str(zero)) )

        return True if shortcircuit else violations

    def print_CoM(self, propagators, formfile):
        for i, vertex in enumerate(self.vertices):
            print(f"local [{self.name}:{vertex.name or i+1}] = {'+'.join(str(p) for p in vertex.ingoing_momenta)};", file=formfile)

    def print_graph(self, propagators, mathfile):
        VertexShapeMap = { 2 : '"Circle"', 4 : '"Circle"', 6 : '"Square"', 8 : '"Triangle"' }
        VertexSizeMap  = { 2 : 0.001, 4 : 0.1, 6 : 0.15, 8 : 0.2 }

        print(indent(dedent(f"""
            {self.name} -> EdgeTaggedGraph[
                {{
                    {newline(5, prefix=',').join(
                        [f"DirectedEdge[{edge.source+1}, {edge.destination+1}, {edge.propagator.momentum}]"
                            for edge in self.edges]
                        +
                        [f"DirectedEdge[mu{ext+1}, {vertex+1}, p{ext+1}]"
                            for ext, vertex in enumerate(self.external)]
                    )}
                }},
                VertexShapeFunction -> {{
                    {newline(5, prefix=',').join(
                        [f'{i+1} -> {VertexShapeMap[vertex.order]}'
                            for i,vertex in enumerate(self.vertices)]
                        +
                        [f'mu{ext+1} -> "Circle"'
                            for ext in range(len(self.external))]
                    )}
                }},
                VertexSize -> {{
                    {newline(5, prefix=',').join(
                        [f'{i+1} -> {VertexSizeMap[vertex.order]}'
                            for i,vertex in enumerate(self.vertices)]
                        +
                        [f'mu{ext+1} -> 0.001'
                            for ext in range(len(self.external))]
                    )}
                }},
                VertexLabels -> None,
                VertexStyle -> Black,
                EdgeStyle -> Black,
                EdgeLabels -> "EdgeTag"
            ]"""), " "*4), file=mathfile, end='')

    def list_vertices(self):
        vertex_roster = {}

        for vertex in self.vertices:
            vertex_roster[vertex] = (vertex_roster.get(vertex, 0)) + 1
        return vertex_roster

    @staticmethod
    def symmetry_factor_FORM(factor):
        return re.sub(r'([0-9]+)!', r'fac_(\1)', factor)

    def permute_FORM(self, perm, legs):
        replace = ', '.join(f'p{i+1},p{perm[i]+1}' for i in range(len(perm)))
        for p in PARTICLES:
            if p not in legs:
                continue

            try:
                restricted = perm.on_range(self.legs_up_to(legs, p) - legs[p], self.legs_up_to(legs, p))
            except ValueError:
                raise ChPTError(f"External leg permutation {perm} exchanges {PARTICLES[p]}s and non-{PARTICLES[p]}s")

            match p:
                case 'M':
                    array = 'extflav'
                case 'V'|'A':
                    array = 'external'
                case _:
                    continue

            #logger.debug(f"Permutation {perm} on {p} -> replace {restricted}", extra=logextra)
            replace += ',\n' + ' '*(4 + len('+ replace_(')) + ', '.join(f'{array}[{i+1}],{array}[{restricted[i]+1}]' for i in range(legs[p]))


        return f"replace_({replace})"

    def define_FORM(self, global_vertex_roster, propagators, legs, loops, formdir, prefix = ''):
        vertex_roster = {}
        vertex_tag_map = {}

        # Figure out how many vertices of each kind are needed in the diagram,
        # and map the local vertex tag (as used by the diagram) to the global one
        # (as used by the set of identical vertices of the given kind defined in FORM)
        for index, vertex in enumerate(self.vertices):
            vertex_tag_map[index] = vertex_roster.get(vertex, 0)
            vertex_roster[vertex] = vertex_tag_map[index] + 1

        # Update the global vertex roster
        for vertex, count in vertex_roster.items():
            if (global_vertex_roster.get(vertex, 0)) < count:
                global_vertex_roster[vertex] = count

        diagram = f"{prefix}diagram{self.name}"

        # Write the diagram's flags
        with open(f"{formdir}/flags/on_{diagram}.hf", 'w') as on, open(f"{formdir}/flags/off_{diagram}.hf", 'w') as off:
            print_info_FORM(on,  f"This file enables the flags for {diagram}")
            print_info_FORM(off, f"This file disables the flags for {diagram}")

            for flag,value in self.flags:
                if value is None:
                    print(f'#define {flag}', file=on)
                else:
                    print(f'#define {flag} "{value}"', file=on)
                print(f'#undefine {flag}', file=off)

            logger.log(VERBOSE, f"Wrote output file {on.name}", extra=logextra)
            logger.log(VERBOSE, f"Wrote output file {off.name}", extra=logextra)

        # Define the diagram as the product of its vertices
        with open(f"{formdir}/diagrams/{diagram}.hf", 'w') as formfile:

            print_info_FORM(formfile, f"This file defines and assembles {diagram}.")

            # Compose diagram
            print(dedent(f"""\
                        global {diagram} = {' * '.join(
                            f'i_*vert{v.name_FORM(vertex_tag_map[i])}'
                            for i,v in enumerate(self.vertices))
                        }"""), file=formfile, end='')

            # Divide by symmetry factor
            if len(self.symmetry_factor) > 2 or self.symmetry_factor[None] != '1':
                print(indent(dedent(f"""
                            #ifdef `NOSYMFACT'{''.join(f'''
                            #elseif isdefined({var})
                                / ({self.symmetry_factor_FORM(sym)})'''
                                for var,sym in self.symmetry_factor.items() if var)}
                            #else
                                / ({self.symmetry_factor_FORM(self.symmetry_factor[None])})
                            #endif
                            ;"""),' '*4), file=formfile)
            else:
                print(';', file=formfile)


            # Connect external lines
            ext = 0
            for particle in PARTICLES:
                for _ in range( legs.get(particle, 0) ):
                    tag = self.vertices[self.external[ext]].name_FORM(vertex_tag_map[self.external[ext]])
                    ext += 1
                    match particle:
                        case 'M':
                            print(f"id,all phi(flav?f{tag}x, ?lorentz) = replace_(flav, extflav[{ext}]) * derivs(p{ext}, ?lorentz);", file=formfile)
                        case 'V':
                            # NOTE: the use of ext is to ensure that Schoonschip notation doesn't spoil the
                            #       pattern matching of subsequent id's, which it would with a naked vector
                            # NOTE: QED uses A for the photon, otherwise that is the axial vector
                            print('\n'.join([
                                "#ifdef `QED'",
                                f"    id,all A(?lorentz, mu?mu{tag}x) = ext(external[{ext}], mu) * derivs(p{ext}, ?lorentz);",
                                "#else",
                                f"    id,all V(?lorentz, mu?mu{tag}x) = ext(external[{ext}], mu) * derivs(p{ext}, ?lorentz);",
                                "#endif"
                                ]), file=formfile)
                        case 'A':
                            print(f"id,all A(?lorentz, mu?mu{tag}x) = ext(external[{ext}], mu) * derivs(p{ext}, ?lorentz);", file=formfile)
                        case 'S' | 'P':
                            print(f"id,all {particle}(?lorentz) = derivs(p{ext}, ?lorentz);", file=formfile)
                        case _:
                            raise ChPTError(f"Particle type '{particle}' not implemented for pickout")

            # Connect the edges, picking out fields from the respective vertices
            for i,edge in enumerate(self.edges):
                src_tag = self.vertices[edge.source     ].name_FORM(vertex_tag_map[edge.source     ])
                dst_tag = self.vertices[edge.destination].name_FORM(vertex_tag_map[edge.destination])
                momentum = edge.propagator.momentum

                if edge.propagator.flav_dependent:
                    idxa = f'index{i+1}a'
                    idxb = f'index{i+1}b'
                    prop = f'propmatrix({edge.propagator.momentum}, {idxa}, {idxb})'
                else:
                    idxa = f'index{i+1}'
                    idxb = idxa
                    prop = f'prop({edge.propagator.momentum},{edge.propagator.mass_squared})'

                print(dedent(f"""\
                    multiply i_ * {prop};
                    id,all phi(flav?f{src_tag}x, ?lorentz) = replace_(flav, {idxa}) * derivs(-({momentum}), ?lorentz);
                    id,all phi(flav?f{dst_tag}x, ?lorentz) = replace_(flav, {idxb}) * derivs(+({momentum}), ?lorentz);"""
                    ), file=formfile)

            # Do all contractions, etc.
            print(dedent(f"""\
                id ext(ext1?, mu?) = ext1(mu);
                #call doderivs
                #call dotrace({diagram})"""), file=formfile)

            logger.log(VERBOSE, f"Wrote output file {formfile.name}", extra=logextra)

        # If necessary, permute external legs
        with open(f"{formdir}/permute/{diagram}.hf", 'w') as formfile:
            if len(self.permutation_group) > 1:
                logger.debug(f"Permutation group of {self.name} [size {len(self.permutation_group)}] is  {', '.join(perm.oneline_string(sep=',', base=1) for perm in self.permutation_group)}", extra=logextra)

                print_info_FORM(formfile, f"This file permutes the external legs of {diagram}.")
                print(dedent(f"""\
                    #call nskip({diagram})
                    multiply (1"""), file=formfile)
                for perm in self.permutation_group:
                    if perm.is_identity():
                        continue
                    print(f"    + {self.permute_FORM(perm, legs)}", file=formfile)
                print(dedent(f"""\
                                );
                            .sort:>>permute {diagram}<<;"""), file=formfile)
            else:
                print_info_FORM(formfile, f"This file is an empty placeholder, since the permutation group is trivial.")

            logger.log(VERBOSE, f"Wrote output file {formfile.name}", extra=logextra)


        with open(f"{formdir}/loops/{diagram}.hf", 'w') as formfile:

            print_info_FORM(formfile, f"This file identifies the loop integrals in {diagram}.")

            # Identify loop integrals FIXME
            loops_present = self.get_loop_momenta(loops)
            if loops_present:
                # Function representing integral over all present loop momenta
                integral = f"i_^{len(loops_present)} * int{''.join(sorted(loops_present))}"
                # All propagators that contain only the present loop momenta
                props_present = [i for i,p in enumerate(propagators)
                                 if (p.momentum.vectors() & loops) <= loops_present]

                print(dedent(f"""\
                            #call nskip({diagram})
                            id {'*'.join(f'prop{i+1}^n{i+1}?' for i in props_present)}
                                = {integral}({','.join(f'n{i+1}' for i in props_present)});
                            """),
                        file=formfile)

            logger.log(VERBOSE, f"Wrote output file {formfile.name}", extra=logextra)

        return diagram

    def get_loop_momenta(self, loops):
        return {l for edge in self.edges for l in edge.propagator.momentum.vectors() if l in loops}

    def get_order(self, loops):
        return 2 + 2*len(self.get_loop_momenta(loops)) + sum(vert.order - 2 for vert in self.vertices)

    def get_particle_content(self):
        return {p for vert in self.vertices for p in vert.legs.keys()}

def specify_legs(token, legs):
    try:
        if token.isnumeric():
            particle = 'M'
            number = int(token)
        elif token[0].isnumeric():
            particle = token[-1]
            number = int(token[:-1])
        elif token in PARTICLES:
            particle = token
            number = 1
        else:
            return False
    except ValueError as err:
        raise ChPTError(f"Failed to read number of {PARTICLES[particle]}s: {err}")

    if particle not in PARTICLES:
        raise ChPTError(f"Invalid particle type: {particle}")
    if particle in legs:
        raise ChPTError(f"Number of {PARTICLES[particle]} legs already specified")
    legs[particle] = number

    return True

class Momentum:

    def __init__(self, momentum, vectors=set(), implicit=False):
        if isinstance(momentum, Momentum):
            self.components = copy(momentum.components)
        elif isinstance(momentum, str):
            self.components = Momentum.parse_momentum(momentum, vectors, implicit)
        elif isinstance(momentum, Mapping):
            self.components = momentum
        elif momentum == 0:
            self.components = {}

    @staticmethod
    def parse_momentum(string, vectors, implicit):
        # This could be done using an auxiliary parser class, but the amount of state is
        #  so puny (just the index and copies of the arguments of this method) that it's
        #  easier to just put the class methods as inner methods and make the necessary
        #  variables nonlocal.

        index = 0

        # Location in string, for error reporting
        def where():
            nonlocal string, index
            return f"in '{string}' (character {index})"

        # Parses an expression up to end-of-string or close-paren
        def get_expression(depth):
            nonlocal string, index

            result = None
            term = get_term(depth)

            # Obtain terms until none remain, add them together
            while term is not None:
                if result is None:
                    result = term
                else:
                    if isinstance(result, Momentum) != isinstance(term, Momentum):
                        #print(f"{term=}")
                        raise ChPTError(f"Attempting to add scalar and vector {where()}")
                    result = result + term

                term = get_term(depth)

            if result is None:
                raise ChPTError(f"Empty expression {where()}")
            if depth == 0 and index < len(string):
                raise ChPTError(f"Unexpected ')' {where()}")

            index += 1
            return result

        # Parses a term (product of factors preceded by a string of +/-
        # The +/- string is always optional but the previous term will not end correctly
        #  unless at least one +/- is used as a separator.
        # Returns None if the expression ends
        def get_term(depth):
            nonlocal string, index, vectors, implicit

            term = Fraction(1)
            has_sign = False
            has_factor = False

            while index < len(string):
                char = string[index]

                if char.isspace():
                    pass
                elif char == '+':
                    has_sign = True
                elif char == '-':
                    has_sign = True
                    term *= -1
                else:
                    # Obtain factors until none remain, multiply them together
                    factor, oper = get_factor(depth, first=True)
                    while factor is not None:

                        has_factor = True

                        if isinstance(factor, Momentum):
                            if oper == '/':
                                raise ChPTError(f"Division by vector {where()}")
                            if isinstance(term, Momentum):
                                raise ChPTError(f"Product of vectors {where()}")
                            term = term * factor
                        elif isinstance(factor, str):
                            if oper == '/':
                                raise ChPTError(f"Division by vector {where()}")
                            if isinstance(term, Momentum):
                                raise ChPTError(f"Product of vectors {where()}")
                            if factor not in vectors:
                                if not implicit:
                                    raise ChPTError(f"Unknown vector: {factor}")
                                vectors.add(factor)
                            term = Momentum({factor: term})
                        else:
                            term = term * (Fraction(factor,1) if (oper == '*') else Fraction(1,factor))

                        factor, oper = get_factor(depth)
                    break

                index += 1

            if not has_factor:
                if has_sign:
                    raise ChPTError(f"Expected term following '+' or '-' {where()}")
                return None

            return term

        # Parses a factor (a number, a momentum, or a parenthesized expression)
        #  and returns it along with the preceding operator (* or /)
        # The operator is an implicit * if the factor is the first in the term
        # Returns None,None if the term (and possibly the expression) ends
        def get_factor(depth, first=False):
            nonlocal string, index

            oper = '*' if first else None

            while True:
                char = string[index] if index < len(string) else 'end-of-string'

                if char.isspace():
                    pass
                elif char == '*' or char == '/':
                    if oper:
                        raise ChPTError(f"Unexpected {char} {where()}")
                    oper = char
                elif char == '(':
                    index += 1
                    return get_expression(depth+1), oper if oper else '*'
                elif char == '+' or char == '-' or char == ')' or char == 'end-of-string':
                    if oper and not first:
                        raise ChPTError(f"Expected factor between '{oper}' and '{char}' {where()}")
                    return None, None
                elif char.isalnum():
                    start = index
                    index += 1
                    while index < len(string) and string[index].isalnum():
                        index += 1

                    # try-except is inefficient, but the only foolproof way I know of figuring out
                    # if a string is a valid int is to try to convert it to one
                    try:
                        factor = int(string[start:index])
                        #print(f"Numeric factor of {factor}")
                    except ValueError:
                        factor = string[start:index].strip()
                        #print(f"Vector factor of {factor}")

                    return factor, oper if oper else '*'
                else:
                    raise ChPTError(f"Invalid operator '{char}' {where()}")

                index += 1

        # End of inner methods.
        # This just wraps the initial call to get_expression

        result = get_expression(depth=0)
        if not isinstance(result, Momentum):
            if result == 0:
                result = Momentum(result)
            else:
                raise ChPTError(f"Expected momentum, but expression is scalar: '{string}'")

        return result.components

    def vectors(self):
        return {vec for vec,coeff in self.components.items() if coeff != 0}

    def substitute(self, substitutions):
        for tries in range(len(substitutions)+1):
            done = True
            for vec,replacement in substitutions:
                if vec in self.components and self.components[vec] != 0:
                    done = False
                    coeff = self.components.pop(vec)
                    self += coeff * Momentum(replacement)
            if done:
                return self
        raise ChPTError(f"Infinite loop of substitutions encountered")
    def substituted(self, substitutions):
        return Momentum(self).substitute(substitutions)

    def __pos__(self):
        return Momentum({vec: +(self.components[vec]) for vec in self.vectors()})
    def __neg__(self):
        return Momentum({vec: -(self.components[vec]) for vec in self.vectors()})
    def __iadd__(self, othr):
        if isinstance(othr, Momentum):
            for vec, coeff in othr.components.items():
                self.components[vec] = self.components.get(vec, Fraction(0)) + coeff
        elif othr != 0:
            raise TypeError(f"Attempting to add scalar to momentum")
        return self
    def __add__(self, othr):
        sum = Momentum(self)
        sum += othr
        return sum
    def __radd__(self, othr):
        return self + othr
    def __isub__(self, othr):
        if isinstance(othr, Momentum):
            for vec, coeff in othr.components().items():
                self.components[vec] = self.components.get(vec, Fraction(0)) - coeff
        elif othr != 0:
            raise TypeError(f"Attempting to add scalar to momentum")
        return self
    def __sub__(self, othr):
        sum = Momentum(self)
        sum -= othr
        return sum
    def __rsub__(self, othr):
        return (self - othr) * -1
    def __imul__(self, othr):
        if isinstance(othr, Momentum):
            raise TypeError(f"Momentum does not allow scalar products")
        self.components = {vec : coeff*othr for vec,coeff in self.components.items()}
        return self
    def __mul__(self, othr):
        prod = Momentum(self)
        prod *= othr
        return prod
    def __rmul__(self, othr):
        return self * othr
    def __idiv__(self, othr):
        if isinstance(othr, Momentum):
            raise TypeError(f"Momentum does not allow division")
        self.components = {vec : coeff*othr for vec,coeff in self.components.items()}
        return self
    def __div__(self, othr):
        quot = Momentum(self)
        quot /= othr
        return quot

    def __eq__(self, othr):
        if isinstance(othr, Momentum):
            svecs = self.vectors()
            ovecs = othr.vectors()
            return svecs == ovecs and all(self.components(v) == othr.components(v) for v in svecs)
        else:
            return othr == 0 and all(coeff == 0 for coeff in self.components.values())

    def __str__(self):
        return (lambda s : (s[1:] if s[0]=='+' else s) if s else '0')(''.join(f"{'-' if coeff < 0 else '+'}{f'{abs(coeff)}*' if abs(coeff) != 1 else ''}{vec}" for vec,coeff in self.components.items() if coeff != 0))
    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(tuple(sorted((vec,self.components[vec]) for vec in self.vectors())))



class Vertex:
    def __init__(self, tokens):
        if not tokens:
            raise ChPTError("Expected vertex definition after 'V'")

        self.order = self.name = None
        self.legs = {}
        self.ingoing_momenta = []

        for token in tokens:
            if specify_legs(token, self.legs):
                pass

            elif token[-2:] == 'LO':
                if self.order:
                    raise ChPTError("Order already specified")
                self.order = order_NtoO(token)

            else:
                if self.name:
                    raise ChPTError("Vertex name already specified")
                if token.strip().isnumeric():
                    raise ChPTError("Vertex name could be confused for vertex index")
                self.name = token

        if not self.order:
            self.order = 2
        if not self.legs:
            raise ChPTError("No number of legs specified for vertex")

    def n_legs(self):
        return sum(self.legs.values())

    def __str__(self):
        return f"{order_OtoN(self.order)}.{'.'.join(f'{n}{p}' for p,n in self.legs.items())}"
    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        # Note that names are not compared
        return self.order == other.order and self.legs == other.legs

    def __hash__(self):
        return hash((self.order, *self.legs.items()))

    def define_FORM(self, formfile, tag):
        if set(self.legs.keys()) > set(PARTICLES.keys()):
            raise ChPTError(f"Particles other than {and_join(list(PARTICLES.values()))} not supported by ChPT_vertex.hf")

        name = self.name_FORM(tag)

        print(f"#call makevertex({name},{self.order},{','.join(str(self.legs.get(p, 0)) for p in PARTICLES)})", file=formfile)

        return f"vert{name}"

    def name_FORM(self, tag=None):
        return f"{order_OtoN(self.order)}x{''.join(f'{p}{n}' for p,n in self.legs.items() if n or p == 'M')}{'' if tag is None else f'x{tag+1}'}"

class Propagator:
    def __init__(self, tokens, momenta):
        if len(tokens) != 2:
            raise ChPTError("Invalid propagator specification: expected momentum and mass")

        self.momentum = Momentum(tokens[0], vectors=momenta, implicit=False)
        self.mass_squared = tokens[1]
        self.flav_dependent = (self.mass_squared == '*')

    def __str__(self):
        return f"({self.momentum})^2 - {'(flav-dependent)' if self.flav_dependent else self.mass_squared}"

class Edge:
    def __init__(self, vertices, vertex_name_map, propagators, propagator_name_map, tokens):
        if len(tokens) != 3 or any(not token for token in tokens):
            raise ChPTError("Invalid edge specification: expected source, destination and propagator indices/names")

        def get_index(token, list, name_map, type):
            if token in name_map:
                index = name_map(token)
            elif token.isnumeric():
                index = int(token) - 1
                if index < 0 or index >= len(list):
                    raise ChPTError(f"Index for {type} out of bounds: {token}")
            else:
                raise ChPTError(f"Unknown {type}: {token}")

            return index

        self.source      = get_index(tokens[0], vertices, vertex_name_map, "vertex")
        self.destination = get_index(tokens[1], vertices, vertex_name_map, "vertex")
        self.unresolved = (tokens[2][0] == '*')
        if self.unresolved:
            self.propagator = tokens[2][1:]
        else:
            self.propagator = propagators[get_index(tokens[2], propagators, propagator_name_map, "propagator")]

    def __str__(self):
        return f"{self.source+1}-[{'(unresolved)' if self.unresolved else self.propagator}]->{self.destination+1}"

DEFAULT_STEPS = {
    'includes': [
        "* [includes]",
        "#include- {0}/definitions.hf",
        "#include- {0}/kinematics.hf"],
    'replarg': [
        "* [replarg] This flag enables replacements also inside function arguments",
        "*     Remove if e.g. propagators are to be handled separately",
        "#define REPLARG"],
    'vertices': [
        "* [vertices] This creates all the vertices, and may be a bit time-consuming.",
        "* The extraction from the Lagrangian is cached, so it should run faster in subsequent runs.",
        "#include- {0}/vertices.hf"],
    'diagram': [
        "* [diagram] This sets up the Feynman rules of the diagram and contracts the flavor indices",
        "    #include- {0}/diagrams/`DIAGRAM'.hf"],
    'permute': [
        "* [permute] This performs all permutations of the diagram's external legs",
        "    #include- {0}/permute/`DIAGRAM'.hf"],
    'replace': [
        "* [replace] This enacts all replacements specified with 'R'",
        "    #call replacements(`DIAGRAM')"],
    'loops': [
        "* [loops] This identifies loop integrals",
        "    #call fullkinematics(`DIAGRAM')",
        "    #include- {0}/loops/`DIAGRAM'.hf"],
    'check': [
        "* [check] Check that all propagators were substituted",
        "    #call checkprops"],
    'postprocess': ["* [postprocess]"],
    'print': [
        "* [print] Print the results",
        "bracket `LECS', F, i_;",
        "print +s;"]
    }
DEFAULT_PROCEDURES = {
    'nskip': ([], []),
    'trivialkinematics': (
        [".sort", "#call nskip(`DIAGRAM')"],
        [".sort:>>trivial kinematics<<;"]
        ),
    'fullkinematics': (
        [".sort", "#call trivialkinematics(`DIAGRAM')", "#call nskip(`DIAGRAM')"],
        [".sort:>>full kinematics<<;", "#call replacements(`DIAGRAM')"]
        ),
    'replacements': (
        [".sort", "#call nskip(`DIAGRAM')"],
        [".sort:>>replacements<<;"]
        ),
    'replacementlist': ("", "")
    }
class Customization:
    def __init__(self):
        self.before = []
        self.after = []
        self.keep_default = True

class DiagramSet:

    def __init__(self, filename, options=[]):

        self.name = None
        self.order = None
        self.diagrams = {}
        self.propagators = []
        self.propagator_name_map = {}
        self.vector_replacements = []
        self.other_replacements = []
        self.scalar_products = {}
        self.loop_momenta = []
        self.independent_momenta = []
        self.independent_products = []
        self.legs = {}
        self.options = options

        self.customizations = {key: Customization() for key in chain(DEFAULT_STEPS, DEFAULT_PROCEDURES)}

        if not Path(filename).is_file():
            filename = f"{filename}.chpt"

        current_diagram = None
        with open(filename, 'r') as file:
            for number, line in enumerate(file):
                current_diagram = self.parse_statement(f"{filename}:{number+1}", line, current_diagram)
        self.add_diagram(current_diagram)
        logextra['where'] = filename

        if not self.name:
            self.name = Path(filename).stem
            logger.warning(f"Unnamed diagram set; 'N {self.name}' implied by filename", extra=logextra)
        if not self.legs:
            raise ChPTError("External legs not specified, please use 'X ...'")
        if not self.order:
            logger.warning("No diagrams specified", extra=logextra)
        if not self.independent_momenta:
            raise ChPTError(f"Independent momenta not specified, suggestion: M {' '.join(f'p{i+1}' for i in range(len(self.legs)-1))}")
        if not self.vector_replacements:
            logger.warning("Momentum replacements not specified, are you sure conservation of momentum is ensured?", extra=logextra)
        if not self.scalar_products:
            logger.warning("Propagator basis relations not specified, can be determined with command 'basis'.", extra=logextra)

        logger.info(f"Read input file {filename}", extra=logextra)

        self.particle_content = {f'{p}({PARTICLES[p][1:]})' for p in set(self.legs) | {p for diagr in self.diagrams.values() for p in diagr.get_particle_content()}}
        logger.info(f"Particle content is {and_join(self.particle_content)}", extra=logextra)


    def parse_statement(self, where, line, diagram):
        logextra['where'] = where

        # Discard commented-out portion of line (overridden in C-statements)
        if not line.startswith('C') and (comment := line.find('#')) != -1:
            line = line[:comment]

        # Tokenize, ignore empty lines
        tokens = [token.strip() for token in re.split(r'\s+', line) if token.strip()]
        if not tokens:
            return diagram

        # Process tokens, dress exceptions with location in file
        try:

            # Enforce preample
            if tokens[0] in 'NXMRLBPC' and diagram is not None:
                raise ChPTError(f"'{tokens[0]}' may only be used in the preamble (before the first 'D')")
            if tokens[0] in 'VESGF' and diagram is None:
                raise ChPTError(f"'{tokens[0]}' may only be used inside a diagram (after 'D')")

            match tokens[0]:

                case 'I':
                    for filename in tokens[1:]:
                        logger.debug(f"Importing file '{filename}'...", extra=logextra)
                        with open(filename, 'r') as file:
                            for n,l in enumerate(file):
                                diagram = self.parse_statement(f"{where}>{filename}:{n+1}", l, diagram)

                        logextra['where'] = where
                        logger.info(f"Read input file {filename}", extra=logextra)

                case 'N':
                    if self.name:
                        raise ChPTError("Name already specified")
                    if len(tokens) != 2:
                        raise ChPTError(f"Invalid name: {' '.join(tokens[1:])}")
                    self.name = tokens[1]

                    logger.debug(f"Named set '{self.name}'", extra=logextra)

                case 'X':
                    self.init_external(tokens[1:])

                case 'M':
                    for token in tokens[1:]:
                        self.independent_momenta.append(token)
                        self.independent_products += [(q,token) for q in self.independent_momenta]
                        logger.debug(f"Added independent momentum '{token}'", extra=logextra)

                case 'R':
                    if len(tokens) != 3:
                        raise ChPTError("Expected exactly two tokens (left- and right-hand side) in a replacement")
                    self.add_replacement(tokens[1], tokens[2])
                    logger.debug(f"Added replacement {tokens[1]} -> {tokens[2]}", extra=logextra)

                case 'L':
                    for token in tokens[1:]:
                        self.loop_momenta.append(token)
                        logger.debug(f"Added loop momentum '{token}'", extra=logextra)

                case 'B':
                    expected_len = 3 + len(self.propagators + self.independent_products)
                    #if len(tokens) != expected_len:
                        #raise ChPTError(f"Invalid number of elements in propagator basis relation ({len(tokens)}, expected {expected_len})")

                    self.scalar_products[(tokens[1],tokens[2])] = tokens[3:]
                    logger.debug(f"""Loaded relation {tokens[1]}*{tokens[2]} = {' '.join(
                            (f'- {coeff[1:]}*' if coeff.startswith('-') else f'+ {coeff}*' if coeff != '1' else '')
                            + (f'({self.propagators[i].momentum})^2'
                                if i < len(self.propagators) else
                                '.'.join(str(q) for q in self.independent_products[i - len(self.propagators)]))
                            for i, coeff in enumerate(tokens[3:]) if coeff != '0')
                        }""", extra=logextra)

                case 'P':
                    prop = Propagator(tokens[1:], self.valid_momenta())

                    if prop.momentum in self.propagator_name_map:
                        raise ChPTError(f"More than one propagator with momentum {prop.momentum}")
                    self.propagator_name_map[prop.momentum] = len(self.propagators)

                    self.propagators.append(prop)
                    logger.debug(f"Added propagator {len(self.propagators)}: {prop}", extra=logextra)

                case 'C':
                    if len(tokens) < 2:
                        raise ChPTError(f"Missing key in C-statement")
                    key = tokens[1]
                    self.customize(key, line[line.find(key)+len(key)+1:].rstrip())

                case 'D':
                    if len(tokens) < 2:
                        raise ChPTError(f"Missing diagram name in D-statement")
                    self.add_diagram(diagram)
                    diagram = Diagram(tokens[1], sum(self.legs.values()))
                    if diagram.name in self.diagrams:
                        raise ChPTError(f"More than one diagram with name {diagram.name}")
                    logger.debug(f"Initialized diagram '{diagram.name}'", extra=logextra)

                case 'V':
                    if len(tokens) >= 2 and tokens[1][0] == '@':
                        reference = tokens[1][1:]
                        if reference not in self.diagrams:
                            raise ChPTError(f"Unknown reference to diagram '{reference}'")

                        logger.debug(f"Importing vertices from diagram '{reference}'...", extra=logextra)
                        for vertex in self.diagrams[reference].vertices:
                            diagram.add_vertex(copy(vertex))
                    else:
                        diagram.add_vertex( Vertex(tokens[1:]) )

                case 'E':
                    if len(tokens) >= 2 and tokens[1][0] == '@':
                        reference = tokens[1][1:]
                        if reference not in self.diagrams:
                            raise ChPTError(f"Unknown reference to diagram '{reference}'")

                        logger.debug(f"Importing edges from diagram '{reference}'...", extra=logextra)
                        for edge in self.diagrams[reference].edges:
                            diagram.add_edge(copy(edge))
                    else:
                        diagram.add_edge( Edge(diagram.vertices, diagram.vertex_name_map, self.propagators, self.propagator_name_map, tokens[1:]) )

                case 'S':

                    if len(tokens) < 2:
                        raise ChPTError(f"Missing symmetry factor in S-statement")

                    diagram.set_symmetry_factor(tokens[1:])

                case 'G':
                    if len(tokens) >= 2 and tokens[1][0] == '@':
                        reference = tokens[1][1:]
                        if reference not in self.diagrams:
                            raise ChPTError(f"Unknown reference to diagram '{reference}'")

                        logger.debug(f"Importing symmetry group from diagram '{reference}'...", extra=logextra)
                        diagram.permutation_group = self.diagrams[reference].permutation_group
                    else:
                        diagram.set_permutation_group(tokens[1:])

                case 'F':
                    diagram.add_flags(tokens[1:])

                case _:
                    raise ChPTError(f"Unknown statement '{tokens[0]}'")

        except IndexError:
            raise ChPTError("Missing compulsory argument\n" + indent(f"at {where}{newline(1)}{line}", ' '*4))

        except (ChPTError, FileNotFoundError) as err:
            errstr = str(err)
            if not errstr.endswith('\n'):
                errstr = errstr + '\n'
            raise ChPTError(errstr + indent(f"at {where}{newline(1)}{line}", ' '*4))

        return diagram

    def init_external(self, tokens):

        for token in tokens:
            if not specify_legs(token, self.legs):
                raise ChPTError("Invalid number-of-legs specification")

        logger.debug(f"External legs: {'.'.join(f'{n}{p}' for p,n in self.legs.items())}", extra=logextra)

    def add_diagram(self, diagram):
        if diagram is None:
            return

        if diagram.permutation_group is None:
            diagram.permutation_group = Trivial(diagram.n_legs)
        diagram.connect_momenta(self.propagators, self.legs, self.valid_momenta())
        self.diagrams[diagram.name] = diagram

        order = diagram.get_order(self.loop_momenta)
        if self.order is not None:
            if self.order != order:
                raise ChPTError(f"Inconsistent orders: both {order_OtoN(self.order)} and {order_OtoN(order)} (diagram {diagram.name}) present")
        else:
            self.order = order

        logger.debug(f"Finalized diagram '{diagram.name}'", extra=logextra)

    def add_replacement(self, lhs, rhs):
        try:
            vec = Momentum(rhs, vectors=self.valid_momenta())
            self.vector_replacements.append( (lhs, vec) )
        except ChPTError as err:
            logger.warning(f"Could not parse {lhs} -> {rhs} as a replacement of known vectors", extra=logextra)
            logger.info(f" This replacement will be used only as-is in generated code and not internally for CoM", extra=logextra)
            logger.debug(f" Reason: {err}", extra=logextra)
            self.other_replacements.append( (lhs, rhs) )

    def customize(self, fullkey, line):
        if fullkey == '-d':
            logger.log(VERBOSE, f"Defining FORM preprocessor variable {line.strip()}", extra=logextra)
            self.options.append(line.strip())
            return

        variant = "after"
        if fullkey.endswith(".before"):
            variant = "before"
            key = fullkey[:-len(".before")]
        elif fullkey.endswith(".after"):
            variant = "after"
            key = fullkey[:-len(".after")]
        elif fullkey.endswith("!"):
            variant = "!"
            key = fullkey[:-1]
        else:
            key = fullkey

        try:
            match variant:
                case "!":
                    self.customizations[key].keep_default = False
                case "before":
                    self.customizations[key].before.append(line)
                case "after":
                    self.customizations[key].after.append(line)
        except KeyError:
            raise ChPTError(f"No such step or procedure: {key}")

    def insert_step(self, key, *args):
        custom = self.customizations[key]

        if logger.level <= VERBOSE:
            for line in custom.before:
                logger.log(VERBOSE, f"Inserting before step {key}: {line}", extra=logextra)
            if not custom.keep_default:
                logger.log(VERBOSE, f"Suppressing step {key}", extra=logextra)
            for line in custom.after:
                logger.log(VERBOSE, f"Inserting after step {key}: {line}", extra=logextra)

        return '\n'.join(chain(
            (line.format(*args) for line in custom.before),
            (line.format(*args) for line in DEFAULT_STEPS[key]) if custom.keep_default else [f"* [step {key} disabled]"],
            (line.format(*args) for line in custom.after)
            ))

    def insert_procedure(self, key, body):
        custom = self.customizations[key]

        if not any((custom.keep_default, custom.before, custom.after)):
            logger.log(VERBOSE, f"Omitting procedure {key} altogether", extra=logextra)
            return f"* [procedure {key} disabled]"

        if logger.level <= VERBOSE:
            for line in custom.before:
                logger.log(VERBOSE, f"Inserting before procedure {key}: {line}", extra=logextra)
            if not custom.keep_default:
                logger.log(VERBOSE, f"Suppressing procedure {key}", extra=logextra)
            for line in custom.after:
                logger.log(VERBOSE, f"Inserting after procedure {key}: {line}", extra=logextra)

        head, tail = DEFAULT_PROCEDURES[key]

        return '\n'.join(chain(
            [f"#procedure {key}(DIAGRAM)"],
            head,
            [''],
            custom.before,
            [body if custom.keep_default else f"* [default body of procedure {key} disabled]"],
            custom.after,
            [''],
            tail,
            ["#endprocedure"]
            ))


    def replacements(self):
        return self.vector_replacements + self.other_replacements

    def external_momenta(self):
        return list(f'p{i+1}' for i in range(sum(self.legs.values())))

    def valid_momenta(self):
        return set(self.external_momenta()) | set(self.independent_momenta) | set(self.loop_momenta) | {f'p{i+1}' for i in range(sum(self.legs.values()))}

    def all_propagators(self):
        return set(self.propagators) | {e.propagator for d in self.diagrams.values() for e in d.edges if not e.unresolved}

    def print_definitions_FORM(self, formfile):
        print(dedent(f'''
            autodeclare index index;
            autodeclare cfunction int;
            cfunctions prop{',propmatrix' if any(p.flav_dependent for p in self.all_propagators()) else ''};
            {f"symbols <prop1,n1>,...,<prop{len(self.propagators)},n{len(self.propagators)}>;"
                if self.propagators else '* (no explicit propagators)'};
            {f"symbols {','.join(f'{msq}' for msq in {p.mass_squared for p in self.all_propagators() if not p.flav_dependent})};"
                if any(not p.flav_dependent for p in self.all_propagators()) else '* (no explicit masses)'}
            {f"vectors {','.join(self.external_momenta())};"
                if self.external_momenta() else '* (no external momenta)'}
            {f"vectors {','.join(self.loop_momenta)};"
                if self.loop_momenta else '* (no loops)'}
            {f"vectors {','.join(self.independent_momenta)};"
                if self.independent_momenta else '* (no independent momenta)'}
            '''), file=formfile)

    def print_kinematics_FORM(self, formfile):

        print('\n\n'.join((
            self.insert_procedure('trivialkinematics',
                '\n\n'.join(indent(dedent(f'''\
                    id prop({+prop.momentum.substituted(self.vector_replacements)},{prop.mass_squared}) = prop{i+1};
                    id prop({-prop.momentum.substituted(self.vector_replacements)},{prop.mass_squared}) = prop{i+1};'''), ' '*4)
                    for i,prop in enumerate(self.propagators) if not prop.flav_dependent)
                ),
            self.insert_procedure('fullkinematics',
                '\n'.join(indent(dedent(f'''\
                    id {p}.{q} = {' '.join(
                        (f'- {coeff[1:]}*' if coeff.startswith('-') else f'+ {coeff}*' if coeff != '1' else '')
                        + (f'(1/prop{i+1} + {self.propagators[i].mass_squared})'
                           if i < len(self.propagators) else
                           f'({".".join(str(q) for q in self.independent_products[i - len(self.propagators)])})')
                        for i, coeff in enumerate(coeffs) if coeff != '0')
                    };'''), ' '*4)
                    for (p,q), coeffs in self.scalar_products.items())
                ),
            self.insert_procedure('replacements',
                indent(dedent(f'''\
                    repeat;
                        #call replacementlist(`DIAGRAM')
                        #ifdef `REPLARG'
                            argument;
                                #call replacementlist(`DIAGRAM')
                            endargument;
                        #endif
                    endrepeat;'''),
                ' '*4)),
            self.insert_procedure('replacementlist',
                    '\n'.join(f'    id {lhs} = {rhs};' for lhs,rhs in self.replacements())
                )
            )), file=formfile)

    def sp_symbols(self):
        return {s : s for s in self.external_momenta() + self.loop_momenta + self.independent_momenta}

    def print_header_FORM(self, formfile):
        print("#-", file=formfile)
        print(f"#appendpath {os.getcwd()}", file=formfile)

        Nf = None

        for opt in self.options:
            if opt[0] == '-':
                continue

            eq = opt.find('=')
            if eq >= 0:
                # if opt[:eq] == 'NF':
                #     Nf = int(opt[eq+1:])
                # elif opt[:eq] == 'PAR' and opt[eq+1:] == 'SQRT':
                #     Nf = "2"
                print(f'#define {opt[:eq]} "{opt[eq+1:]}"', file=formfile)
            else:
                print(f"#define {opt}", file=formfile)

        if 'V' in self.particle_content or 'A' in self.particle_content:
            print("#define HASVECTOR", file=formfile)
        if 'S' in self.particle_content or 'P' in self.particle_content:
            print("#define HASSCALAR", file=formfile)

        print(dedent(f'''\
                #define NEXT "{sum(self.legs.values())}"
                #define ORDER "{order_OtoN(self.order)}"
                #define LOOPMOMENTA "{",".join(self.loop_momenta)}"
                #define INTMOMENTA "{",".join(self.independent_momenta)}"
                #include- ChPTlib.hf'''), file=formfile)

    def finalize_diagrams_FORM(self, vertex_roster, formfile):
        print(dedent(f"""\
                .sort:>>all {self.name} diagrams<<;
                off statistics;
                drop {', '.join(f'vert{vertex.name_FORM(tag)}' for vertex,count in vertex_roster.items() for tag in range(count))};
                """), file=formfile)

    def check_CoM(self):
        logger.info("Checking conservation of momentum...", extra=logextra)
        conservation = True
        for diagram in self.diagrams.values()   :
            violations = diagram.check_CoM(self.propagators, self.sp_symbols(), self.replacements())
            for vertex,perm,total in violations:
                logger.error(f"CoM violation in diagram {diagram.name}, vertex {vertex}, permutation {perm}: momenta sum to {total}", extra=logextra)
                conservation = False

        if conservation:
            logger.info("CoM satisfied in all diagrams.", extra=logextra)

    def check_CoM_FORM(self):
        filename = f"{self.name}_CoM.frm"
        with open(filename, 'w') as formfile:
            print("#-\noff statistics;", file=formfile)
            self.print_definitions_FORM(formfile)

            for diagram in self.diagrams.values():
                diagram.print_CoM(self.propagators, formfile)

            for lhs,rhs in self.replacements():
                print(f"id {lhs} = {rhs};", file=formfile)

            print(dedent("""\
                .sort
                #define RESULT "0"
                #do EXPR={`ACTIVEEXPRNAMES_'}
                    #if (termsin(`EXPR') > 0)
                        #write "Diagram:vertex `EXPR' violates CoM by %E",`EXPR'
                        #redefine RESULT "{`RESULT'+1}"
                    #endif
                #enddo
                #if `RESULT' == 0
                    #write "Conservation of momentum is satisfied.\\n"
                #else
                    #write "\\nConservation of momentum failed in `RESULT' places.\\n"
                #endif
                .end
                """), file=formfile)

            logging.info(f"Wrote output file {formfile.name}", extra=logextra)

        logging.info(dedent(f"""\
                     To check conservation of momentum, run the following:
                     {' '*len('[WARNING]')} $ form {filename}"""), extra=logextra)

    def draw_graphs(self):
        filename = f"{self.name}_graphs.m"

        with open(filename, 'w') as mathfile:
            print(dedent(f"""\
                (* {AUTHOR}
                 *
                 * Mathematica file generated by ChPT.py on {datetime.now().strftime("%c")}
                 *  with arguments '{' '.join(sys.argv[1:])}'
                 * Contains a (diagram name) -> EdgeTaggedGraph association, which Mathematica
                 *  should represent as very ugly Feynman diagrams.
                 *)

                {self.name}graphs = <|"""), file=mathfile, end='')
            first = True
            for diagr in self.diagrams.values():
                if not first:
                    print(',', file=mathfile, end='')
                else:
                    first = False
                diagr.print_graph(self.propagators, mathfile)
            print("|>;", file=mathfile)

            logging.info(f"Wrote output file {mathfile.name}", extra=logextra)

        logging.info(f"To view the graphs, load that file into a Mathematica notebook", extra=logextra)

    def determine_prop_basis(self):
        vector_prods = (
            [(self.loop_momenta[i], self.loop_momenta[j]) for i in range(len(self.loop_momenta)) for j in range(i+1)]
            +
            [(l, q) for l in self.loop_momenta for q in self.independent_momenta])

        expressions = (
            [sp.expand(parse_expr(f"({prop.momentum})**2")) for prop in self.propagators]
            +
            [sp.expand(parse_expr(f"{p}*{q}")) for p,q in self.independent_products])

        if len(vector_prods) != len(self.propagators):
            raise ChPTError("Invalid number of propagators in basis")

        try:
            matrix = sp.Matrix([
                [expr.coeff(parse_expr(f"{p}*{q}"))
                    for p,q in vector_prods + self.independent_products]
                for expr in expressions]) ** -1
        except:
            raise ChPTError("Basis of propagators is incomplete or overcomplete")

        basis = [f"B {p} {q} {' '.join(str(matrix[row,col]) for col in range(len(vector_prods + self.independent_products)))}"
                    for row, (p,q) in enumerate(vector_prods)]

        for elem in basis:
            self.parse_statement('', elem, None)

        print(fill(width=80, text=dedent(f"""\
            The following is the projection of scalar products onto the basis of inverse propagators.
            Following each 'B' are the two vectors in the product, then {len(vector_prods)} numbers, then {len(self.independent_products)} more. The ith number among the {len(vector_prods)} first represents the coefficient in front of the squared momentum of propagator i (i.e., the inverse propagator plus the squared mass) in the scalar product. The remaining numbers represent multiples of the products {', '.join(f'{p}.{q}' for p,q in self.independent_products)} in that order. This can be pasted into the input file for future use.""")))
        for elem in basis:
            print(elem)

    def generate_FORM(self):
        dirname = f"{self.name}"
        os.makedirs(f"./{dirname}", exist_ok=True)
        os.makedirs(f"./{dirname}/diagrams", exist_ok=True)
        os.makedirs(f"./{dirname}/loops", exist_ok=True)
        os.makedirs(f"./{dirname}/flags", exist_ok=True)
        os.makedirs(f"./{dirname}/permute", exist_ok=True)

        Nf = None
        vertex_roster = {}

        # Have each diagram write its definitions, accumulating a total roster of vertices needed
        diagram_names = [
            diagram.define_FORM(vertex_roster, self.propagators, self.legs, set(self.loop_momenta), dirname, prefix = f"{self.name}")
            for diagram in self.diagrams.values()
            ]

        # Write the definitions of all vertices to a separate file
        with open(f"{dirname}/vertices.hf", 'w') as formfile:
            print_info_FORM(formfile, "This file defines all vertices needed for the diagrams.")

            vertex_names = [
                vertex.define_FORM(formfile, index)
                for vertex,count in vertex_roster.items()
                for index in range(count)
                ]

            logger.log(VERBOSE, f"Wrote output file {formfile.name}", extra=logextra)

        with open(f"{dirname}/definitions.hf", 'w') as formfile:
            print_info_FORM(formfile, "This file defines all variables specific to this calculation.")

            self.print_definitions_FORM(formfile)

            print(dedent(f'''\
                        * Allow user to override which diagrams are to be defined
                        #ifndef `DIAGRAMNAMES'
                            #define DIAGRAMNAMES "{",".join(diagram_names)}"
                        #endif
                        #define VERTEXNAMES "{",".join(vertex_names)}"
                        '''), file=formfile)

            logger.log(VERBOSE, f"Wrote output file {formfile.name}", extra=logextra)

        with open(f"{dirname}/kinematics.hf", 'w') as formfile:
            print_info_FORM(formfile, "This file defines procedures for simplifying the kinematics.")
            self.print_kinematics_FORM(formfile)

            logger.log(VERBOSE, f"Wrote output file {formfile.name}", extra=logextra)

        logger.info(f"Wrote output directory {dirname}/", extra=logextra)
        return dirname

    def generate_FORM_main(self):
        dirname = self.generate_FORM()
        filename = f"{dirname}.frm"

        with open(filename, 'w') as formfile:

            # Creation message
            print_info_FORM(formfile, dedent(f"""\

                This is the main FORM file.
                Run with form -l {filename} (or tform -wNTHREAD to taste).
                Results are written to 'save/{self.name}.sav'.
                """), auto=False)

            # Define options, include relevant headers
            self.print_header_FORM(formfile)

            print('\n'.join((
                f"* This flag lets various procedures know that the program was generated by ChPT.py",
                f"#define GENERATEDBYCHPTPY",
                f"",
                self.insert_step('includes', dirname),
                f"",
                self.insert_step('replarg', dirname),
                f"",
                self.insert_procedure('nskip', dedent(f"""\
                    * If diagrams are split into multiple expressions
                    *  (say, projected onto different tensor structures),
                    *  modify this so it nskips all variants of DIAGRAM's name.
                        skip; nskip `DIAGRAM';""")),
                f"",
                self.insert_step('vertices', dirname),
                f"",
                f"* Main loop - this will create and process each diagram in turn",
                f"* Each repetition is done with the active diagram's flags enabled.",
                f"#do DIAGRAM={{`DIAGRAMNAMES'}}",
                f"    .sort",
                f"    skip;",
                f"    #include- {dirname}/flags/on_`DIAGRAM'.hf",
                f"",
                f"*    NOTE: any custom manipulations inserted between either of these steps should",
                f"*    be guarded with #call nskip(`DIAGRAM') at the start of each module",
                f"",
                self.insert_step('diagram', dirname),
                f"",
                self.insert_step('permute', dirname),
                f"",
                self.insert_step('replace', dirname),
                f"",
                self.insert_step('loops', dirname),
                f"",
                self.insert_step('check', dirname),
                f"",
                f"    #include- {dirname}/flags/off_`DIAGRAM'.hf",
                f"#enddo",
                f"",
                f".sort:>>all {self.name} diagrams<<;",
                f"off statistics;",
                f"drop `VERTEXNAMES';",
                f"",
                self.insert_step('postprocess', dirname),
                f"",
                self.insert_step('print', dirname),
                f".sort",
                f"#call save({self.name})",
                f".end"
                )), file=formfile)
            logger.info(f"Wrote output file {formfile.name}", extra=logextra)

        logging.info(dedent(f"""\
                     To compute the diagrams, run the following:
                     {' '*len('[WARNING]')} $ form {filename}"""), extra=logextra)

def print_info_FORM(formfile, extra=[], auto=True):
    print(dedent(f"""\
        *** {AUTHOR}
        *** FORM file generated by ChPT.py on {datetime.now().strftime("%c")}
        ***  with arguments '{' '.join(sys.argv[1:])}'"""), file=formfile)
    if auto:
        print(dedent(f"""\
            *** This file is automatically generated and should preferably not be modified.
            *** Instead, the main FORM file is meant to accommodate customization."""), file=formfile)
    for line in extra.split('\n'):
        print(f"*** {line}", file=formfile)

    print("", file=formfile)

def print_chpt_help():
    def wrap(char = None):
        width = 80
        if char is None:
            return TextWrapper(width=width, initial_indent='', subsequent_indent=' ')
        else:
            return TextWrapper(width=width, initial_indent=f' {char} ', subsequent_indent='    ')
    print('\n'.join([
          wrap().fill("The .chpt file format is a very simple way of describing multiloop Feynman diagrams in ChPT using a formulation suitable for master integral reduction. All non-empty lines contain a single statement, consisting of a single character specifying the type of statement followed by a number of whitespace-separated arguments. The character '#' will cause the remainder of the lineto be ignored (except in C-statements, which are read verbatim)."),
          "The statement types are the following:",
          wrap('I').fill("Input the contents of other .chpt files given as arguments."),
          wrap('N').fill("Name the diagram set according to the single argument. If no N-statement is given, the set is named after the .chpt file. The name should only consist of letters (no numbers, underscores, etc.) in order to be compatible with both FORM and Mathematica variable name formats."),
          wrap('X').fill(f"Specify the number of external legs (all diagrams have the same). Each argument is a number followed by a character specifying the type of particle. The particle types are {and_join([f'{p}({n[1:]})' for p,n in PARTICLES.items()])}. The particle specifier defaults to M(eson) if omitted. The number defaults to 1 if omitted. The number of each type of particle may be specified only once, and if it is not specified at all, it is set to zero.  External legs are assigned incoming momenta pi, with i being a 1-based index, assigned incrementally with mesons first, then vectors, etc."),
          wrap('L').fill("Define loop momenta, each one given as a separate argument. Bear in mind that l1, l2, etc. may clash with the names of the LECs."),
          wrap('M').fill("Define external momenta, each one given as a separate argument. All independent external momenta should be explicitly given this way."),
          wrap('R').fill("Specify a replacement to be applied to the expressions. For example, this may be used to implement conservation of external momentum, or to replace the automatically assigned external momenta with user-defined ones. Takes two arguments, the left- and right-hand side of the replacement, respectively. Example: 'R p4 -(p1+p2+p3)' for 4-point diagrams. Bear in mind that the expressions should work in both FORM and Mathematica."),
          wrap('P').fill("Define a propagator. Takes two arguments: the momentum and mass-squared of the propagator. Propagators are numbered in the order they are specified, starting at 1. Bear in mind that for master integral reduction, all scalar products of loop momenta, and all scalar products of a loop momentum and an independent external momentum, must be expressible as a linear combination of squared propagator momenta. Thus, the number of propagators should be E*L + L*(L+1)/2 with L loops and E independent external momenta; this may require the inclusion of propagators that don't actually appear as a propagator in any diagram. The command 'basis' checks that this has been done correctly. The mass-squared may be given as '*', indicating flavor-dependent mass."),
          wrap('B').fill("Specify how a product of momenta is expressed as a linear combination of squared propagator momenta, as described under 'P'. The first two arguments are the momenta being multiplied together, and the following are the coefficients in front of each squared propagator momentum, in the order the propagators were defined. The command 'basis' automatically generates all necessary B-statements."),
          wrap('C').fill(f"Specify a customization of the generated FORM files. The first argument (the KEY) should be the name of a computation step in the main FORM file ({and_join(list(DEFAULT_STEPS))}), each marked by [KEY] in the generated FORM code, or the name of a procedure ({and_join(list(key for key in DEFAULT_PROCEDURES))}), or -d. The latter case is special, and acts like the -d command line option. In the other cases, KEY may be followed by a suffix: KEY.before, KEY.after, or KEY! (lack of a suffix is equivalent to KEY.after). The remainder of the line is then taken as verbatim FORM code (except that in the steps, the sequence '{{0}}' is replaced by the path to the generated FORM directory) and inserted before or after (as specified) the computation step or procedure body. If KEY! is given, the default implementation of the step or procedure body is omitted. If the default procedure body is omitted and not replaced by a custom one, the entire procedure is omitted and should instead be provided elsewhere, e.g., through a separate KEY.prc file. Any number of customizations may be applied to the same KEY and will be applied in the order they are entered."),
          wrap('D').fill("Initialize a new active diagram, and finalize the previous one, if any. [VESGF]-statements act on an active diagram, so a D-statement must be made before these can be used. It is recommended to give all other types of statements before the first D-statement, as a preamble. The compulsory argument gives the name of the diagram; this is combined with the diagram set name to produce the full diagram name."),
          wrap('V').fill("Define a vertex for the active diagram. There are three types of arguments, and they can be given in any order. If the argument is of the form of zero or more N's followed by 'LO', or alternatively 'N', a number representing the number of N's, and 'LO', that specifies the order of the vertex, defaulting to LO if omitted. If the argument would be accepted by an X-statement, it specifies the number of legs on the vertex as described there. Otherwise, the argument specifies the optional name of the vertex. Alternatively, there may be a single argument, consisting of '@' and the name of a previously defined diagram; this imports all vertices of that diagram as if all its V-statements were pasted in place of the 'V @...' Vertices are numbered in the order they are defined, starting at 1."),
          wrap('E').fill("Define an edge for the active diagram. There are three arguments: the number or name of the source vertex, the number or name of the destination vertex, and the number or momentum of the propagator represented by the edge. The momentum of the propagator flows from the source to  the destination. The propagator may instead be given as '*MASS', in which case the momentum will be deduced from conservation. MASS should be the mass-squared of the propagator, or another '*' to indicate flavor-dependent mass. Alternatively, there may be a single argument, consisting of '@' and the name of a previously defined diagram; this imports all edges of that diagram similarly to 'V @...'. Note that vertex names and numbers may be different in this diagram, so this should be used with caution."),
          wrap('S').fill("Endow the current diagram with a symmetry factor. The first argument should be an expression for the symmetry factor, by which the diagram will be divided. Exponents (^) and factorials (!) may be used and will automatically be translated to FORM and Mathematica syntax. Optional symmetry factors may be given, prefixed by the name of a FORM preprocessor variable and a colon. If that variable is defined at the moment the diagram is created, that factor will be used instead (the first one in the order given here takes precedence). There is an implicit NOSYMFACT:1 that overrides other alternative symmetry factors."),
          wrap('G').fill("Endow the current diagram with a group of permutations of its external legs, over which it should be summed. There can be at most one G-statement per diagram. The default arrangement of external momenta is that all momenta belonging to each particle type are assigned in the order they appeared in the X-statements, and for each type, they are assigned to all available legs on each vertex in the order they appeared in the V-statements. The group can be specified in either of three ways: one or more permutations, a distribution, or as a copy of another diagram's group with 'G @...'. Permutations are given in cycle notation [example: (1,2)(5,6,7,8)], and the group will be all distinct permutations that can be formed by composing the specified permutations (including the identity). A distribution is given as a bracket-enclosed list (example: [2,2,4]) whose elements sum to the number of legs of the diagram. The group will be all ways of distributing the external momenta into sets with the given list of sizes, without regard to the order within each set, and without regard to the relative order of sets of equal size. Each set-size may be optionally suffixed by an alphanumeric label, and in that case the relative order is respected between sets with different labels; the special label '*' is considered different from all labels including itself. Normally, each set will correspond to the number of external legs on each vertex, and the group will be all ways of assigning the external legs, up to symmetries of the diagram which can be accounted for with the labels; for instance, [2,2] is appropriate for a simple lollipop diagram with four identical external legs, two on each vertex; [2a,2b] or [2*,2*] is appropriate if the vertices are not interchangeable."),
          wrap('F').fill("Equip a diagram with one or more flags. By 'flags' is meant FORM preprocessor variables that will be defined only while processing the current diagram, for convenience when controlling any diagram- specific procedures. The syntax is as for FORM -d options, i.e. NAME(=VALUE).")]))

def setup_logging(args):

    if args.quiet:
        level=logging.ERROR
        format="[%(levelname)s] %(message)s"
    elif args.debug or args.verbose:
        level=logging.DEBUG if args.debug else VERBOSE
        format="[%(levelname)7s/%(where)s] %(message)s"
    else:
        level=logging.INFO
        format="[%(levelname)7s] %(message)s"

    handler = logging.StreamHandler()
    if sys.stdout.isatty():
        handler.setFormatter(ColorFormatter(format))
    else:
        handler.setFormatter(logging.Formatter(format))

    logging.basicConfig(level=level, handlers=[handler])

def main():

    parser = argparse.ArgumentParser(
        prog = 'ChPT',
        formatter_class=argparse.RawTextHelpFormatter,
        usage="ChPT [-hH] [-vDq] [-cCbgfF]  [-d NAME[=VALUE]] ... FILE",
        description=dedent(f"""\
            {AUTHOR} (version 0.3)
            FORM code generator for chiral perturbation theory diagrams."""))

    parser.add_argument('-H', '--chpt-help', action='store_true',
                        help="print information about the .chpt file format and exit")
    parser.add_argument('-c', '--check-CoM', action='append_const', const='c', dest='actions',
                        help="check that conservation of momentum (CoM) is respected by all input diagrams")
    parser.add_argument('-C', '--check-CoM-form', action='append_const', const='C', dest='actions',
                        help=f"like -c, but generate a FORM program to do the checking{newline()}This is useful if CoM relies on additional processing")
    parser.add_argument('-b', '--compute-basis', action='append_const', const='b', dest='actions',
                        help=f"determine how any scalar product of loop and external momenta can be expressed in terms of the basis of inverse propagators provided{newline()}This is both loaded into the program and printed as a set of B statements, which can be pasted into the .chpt file for future use")
    parser.add_argument('-g', '--generate-graphs', action='append_const', const='g', dest='actions',
                        help="produce Mathematica code for drawing (very ugly) graphs of the diagrams, for verifying that the topologies have been properly implemented")
    parser.add_argument('-f', '--generate-form', action='append_const', const='f', dest='actions',
                        help="produce FORM code for computing the diagrams")
    parser.add_argument('-F', '--generate-form-main', action='append_const', const='F', dest='actions',
                        help=f"produce a customizable FORM template program that uses the files generated with --generate-form{newline()}This implies --generate-form")

    parser.add_argument('-d', '--define', action='append',
                        help=f"define FORM preprocessor variables and pass them to FORM.{newline()}Uses the same -d NAME(=VALUE) syntax as FORM")

    parser.add_argument('-v', '--verbose', action='store_true',
                        help="say more about what is being done, including place of origin for all printouts")
    parser.add_argument('-D', '--debug', action='store_true',
                        help="say even more about what is being done; implies --verbose")
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="suppress all output except what is critical.")

    parser.add_argument('file', type=str, metavar='FILE', nargs='?', default=None,
                        help="the main .chpt input file.")

    args = parser.parse_args()

    # Cheat to endow -H with the same behavior as -h
    if args.chpt_help:
        print_chpt_help()
        return 0
    elif not args.file:
        print(parser.usage, file=sys.stderr)
        print("Error: the following arguments are required: FILE")
        return 1

    setup_logging(args)

    try:
        diagrs = DiagramSet(args.file, options=args.define if args.define else [])

        if args.actions is None:
            raise ChPTError(f"No action specified (any of -cCgbfF)")
        for action in args.actions:
            match action:
                case 'c':
                    diagrs.check_CoM()
                case 'C':
                    diagrs.check_CoM_FORM()
                case 'g':
                    diagrs.draw_graphs()
                case 'b':
                    diagrs.determine_prop_basis()
                case 'f':
                    diagrs.generate_FORM()
                case 'F':
                    diagrs.generate_FORM_main()

        return 0

    except FileNotFoundError as err:
        logger.error(f".chpt file not found: '{args.file}'", extra=logextra)
    except ChPTError as err: # User error
        if args.debug:
            logging.error(traceback.format_exc(), extra=logextra)
        else:
            logging.error(err, extra=logextra)
    except Exception: # Debugging error
        logging.error(traceback.format_exc(), extra=logextra)

    return 1

if __name__ == '__main__':
    retval = main()
    sys.exit(retval)
