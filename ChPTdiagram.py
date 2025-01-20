#!/home/mssjo/python/bin/python3

import re
import sys,os
import subprocess
import logging

from collections.abc import Mapping
from copy import copy, deepcopy
from datetime import datetime
from fractions import Fraction
from pathlib import Path
from textwrap import indent, dedent

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

from permute import Permutation, Span, Sn, Trivial, Distributions

logger = logging.getLogger("ChPTdiagram")

PARTICLES = {'M': 'meson', 'V': 'vector'} # TODO: 'A': 'axial-vector', 'S': 'scalar', 'P': 'pseudoscalar'

# Since \n is not allowed in f-strings, use this instead.
# Supports various indentation options as well for convenience.
def newline(indent=0, *, prefix='', infix='', suffix='', indent_str=' '*4):
    return prefix + '\n' + infix + indent_str*indent + suffix

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

class ChPTDiagram:
    def __init__(self, name, n_ext):

        logger.debug(f"Initializing {n_ext}-leg diagram {name}")
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

    def add_vertex(self, vertex, where=''):
        if vertex.name:
            if vertex.name in self.vertex_name_map:
                raise ChPTError(f"Ambiguous vertex name: '{vertex.name}'")
            self.vertex_name_map[vertex.name] = len(self.vertices)

        self.vertices.append(vertex)
        logger.debug(f"{where} Added vertex {len(self.vertices)}{f' ({vertex.name})' if vertex.name else ''}: {vertex}")

        # This should not carry over when copying a vertex,
        # and it's easier to zero it here than to mess with copy()
        vertex.ingoing_momenta = []

    def add_edge(self, edge, where=''):
        self.edges.append(edge)
        logger.debug(f"{where} Added edge {edge}")

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
                    logger.debug(f"Connected external leg {len(self.external)} ({PARTICLES[particle]}) to vertex {i+1}")
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
                    logger.debug(f"Resolved momentum routing for edge {edge}")

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
                logging.debug(f"CoM, diagram {self.name}, vertex {i+1}, permutation {perm}: {zero}")
                zero.substitute({f'p{j+1}' : f'p{perm[j]+1}' for j in range(len(perm))})
                zero.substitute(replacements)
                if zero != 0:
                    if shortcircuit:
                        return False
                    violations.append( (i+1, perm, str(zero)) )

        return True if shortcircuit else violations

    #def print_CoM(self, propagators, formfile):
        #for i, vertex in enumerate(self.vertices):
            #print(f"local [{self.name}:{vertex.name or i+1}] = {''.join(vertex.ingoing_momenta)};", file=formfile)

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

            #logger.debug(f"Permutation {perm} on {p} -> replace {restricted}")
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

        diagram = f"diagram{prefix}{self.name}"

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
                            print(f"id,all A(?lorentz, mu?mu{tag}x) = ext(external[{ext}], mu) * derivs(p{ext}, ?lorentz);", file=formfile)
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

        # If necessary, permute external legs
        with open(f"{formdir}/permute/{diagram}.hf", 'w') as formfile:
            if len(self.permutation_group) > 1:
                logger.debug(f"Permutation group of {self.name} [size {len(self.permutation_group)}] is  {', '.join(perm.oneline_string(sep=',', base=1) for perm in self.permutation_group)}")

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

        return diagram

    def get_loop_momenta(self, loops):
        return {l for edge in self.edges for l in edge.propagator.momentum.vectors() if l in loops}

    def get_order(self, loops):
        return 2 + 2*len(self.get_loop_momenta(loops)) + sum(vert.order - 2 for vert in self.vertices)

def specify_legs(token, legs):
    if token.isnumeric():
        if 'M' in legs:
            raise ChPTError("Number of meson legs already specified")
        legs['M'] = int(token)

    elif token[0].isnumeric():
        particle = token[-1]
        if particle not in PARTICLES:
            raise ChPTError(f"Invalid particle type: {particle}")
        if particle in legs:
            raise ChPTError(f"Number of {PARTICLES[particle]} legs already specified")
        legs[particle] = int(token[:-1])

    else:
        return False

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

    #class MomentumParser:

        #def __init__(self, string, vectors, implicit):
            #self.string = string
            #self.vectors = vectors
            #self.implicit = implicit

            ##print(f"{self.vectors=} ({self.implicit=})")

            #self.index = 0


        #def where(self):
            #return f"in '{self.string}' (character {self.index})"





        #def _parse(self, depth):

            ##print(f"Parsing '{self.string}' ({depth=})")

            #while self.index < len(self.string):
                #char = self.string[self.index]

                ##print(f"{depth=} {self.index=}[{char}] {self.tmark=}[{self.string[self.tmark]}] {self.fmark=}[{self.string[self.fmark]}] {self.term=} {self.oper=}")

                #if char.isspace():
                    #if self.tmark == self.index:
                        #self.tmark += 1
                    #if self.fmark == self.index:
                        #self.fmark += 1
                #elif char == ')':
                    #if depth == 0:
                        #raise ChPTError(f"Unexpected ')' in momentum '{self.string}'")
                    #break
                #elif char == '(':
                    #self.conclude_factor(may_be_empty=True)
                    #self.index += 1
                    #idx,par = type(self)(self.string[self.index:], self.vectors, self.implicit)._parse(depth+1)
                    #self.incorporate_factor(par)
                    #self.index += idx
                    #self.fmark = self.index+1
                    #self.oper = '*'
                #elif char == '+' or char == '-':
                    #self.conclude_term(may_be_empty=True)
                    #self.oper = '*'
                    #self.tmark = self.index+1
                    #self.fmark = self.index+1
                    #if char == '-':
                        #self.term *= -1
                #elif char == '*' or char == '/':
                    #self.conclude_factor(may_be_empty=False)
                    #self.oper = char
                    #self.fmark = self.index+1

                #self.index += 1

            #self.conclude_term(may_be_empty=False)

            #if depth != 0 and self.index == len(self.string):
                #raise ChPTError(f"Unmatched '(' in momentum '{self.string}'")
            #if depth == 0 and self.index != len(self.string):
                #raise ChPTError(f"Unexpected ')' in momentum '{self.string}'")

            ##print(f"Result of parsing: {self.result}")
            #return self.index, self.result

        #def incorporate_factor(self, factor):
            ##print(f"{self.vectors=}")
            #if isinstance(factor, Momentum):
                #if self.oper == '/':
                    #raise ChPTError(f"Division by vector in {self.string}")
                #if isinstance(self.term, Momentum):
                    #raise ChPTError(f"Product of vectors in momentum {self.string}")
                #self.term = self.term * factor
            #elif isinstance(factor, str):
                #if self.oper == '/':
                    #raise ChPTError(f"Division by vector in {self.string}")
                #if isinstance(self.term, Momentum):
                    #raise ChPTError(f"Product of vectors in momentum {self.string}")
                #if factor not in self.vectors:
                    #if not self.implicit:
                        #raise ChPTError(f"Unknown vector: {factor}")
                    #self.vectors.add(factor)
                #self.term = Momentum({factor: self.term})
            #else:
                #self.term = self.term * (Fraction(factor,1) if (self.oper == '*') else Fraction(1,factor))

        #def conclude_factor(self, may_be_empty=False):

            ## Empty factors are only OK in empty terms (implicit 1)
            #if self.fmark == self.index:
                #if may_be_empty:
                    ##print(f"Empty factor")
                    #return
                #else:
                    #raise ChPTError(f"Number or vector expected in {self.string[self.fmark:self.index]}")

            ## try-except is wasteful, but the only foolproof way I know of figuring out
            ## if a string is a valid number is to try to convert it
            #try:
                #fac = int(self.string[self.fmark:self.index])
                ##print(f"Numeric factor of {fac}")
                #self.incorporate_factor(fac)

            #except ValueError:
                #fac = self.string[self.fmark:self.index].strip()
                ##print(f"Vector factor of {fac}")
                #self.incorporate_factor(fac)


        #def conclude_term(self, may_be_empty=False):

            #if self.index == self.tmark:
                #if may_be_empty:
                    #return
                #else:
                    #raise ChPTError(f"Empty term in {self.string}")

            #self.conclude_factor(may_be_empty=True)
            #if isinstance(self.result, Momentum) and not isinstance(self.term, Momentum):
                ##print(f"{self.term=}")
                #raise ChPTError(f"Missing vector in term: {self.string[self.tmark:self.index]}")

            #self.result = self.result + self.term
            #self.term = Fraction(1)

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
        if set(self.legs.keys()) > {'M', 'V'}:
            raise ChPTError("Particles other than mesons and vectors not supported by ChPTdiagram_vertex.hf")

        name = self.name_FORM(tag)

        print(f"#call makevertex({name},{self.legs.get('M', 0)},{self.legs.get('V', 0)},{self.order})", file=formfile)

        return f"vert{name}"

    def name_FORM(self, tag=None):
        return f"P{self.order}{''.join(f'{p}{n}' for p,n in self.legs.items())}{'' if tag is None else f'x{tag+1}'}"

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

class ChPTDiagramSet:

    def __init__(self, filename):

        self.name = None
        self.order = None
        self.diagrams = {}
        self.propagators = []
        self.propagator_name_map = {}
        self.replacements = []
        self.scalar_products = {}
        self.loop_momenta = []
        self.independent_momenta = []
        self.independent_products = []
        self.legs = {}

        if not Path(filename).is_file():
            filename = f"{filename}.xpt"

        current_diagram = None
        with open(filename, 'r') as file:
            for number, line in enumerate(file):
                current_diagram = self.parse_statement(f"{filename}:{number+1}", line, current_diagram)
        self.add_diagram(current_diagram)

        if not self.name:
            raise ChPTError("Diagram set name not specified, please use 'N NAME'")
        if not self.legs:
            raise ChPTError("External legs not specified, please use 'X ...'")
        if not self.order:
            raise ChPTError("No diagrams specified, plese use 'D ...'")
        if not self.independent_momenta:
            raise ChPTError(f"Independent momenta not specified, suggestion: M {' '.join(f'p{i+1}' for i in range(len(self.legs)-1))}")
        if not self.replacements:
            logger.warning("Momentum replacements not specified, are you sure conservation of momentum is ensured?")
        if not self.scalar_products:
            logger.warning("Propagator basis relations not specified, can be determined with command 'basis'.")

        logger.info(f"Read input file {filename}")

    def parse_statement(self, where, line, diagram):
        tokens = [token.strip() for token in re.split(r'\s+', line) if token.strip()]

        # Detect empty lines and full-line comments
        if not tokens or not tokens[0] or tokens[0].startswith('#'):
            return diagram

        try:
            match tokens[0]:

                case 'N':
                    if self.name:
                        raise ChPTError("Name already specified")
                    if len(tokens) != 2:
                        raise ChPTError(f"Invalid name: {' '.join(tokens[1:])}")
                    self.name = tokens[1]

                    logger.debug(f"{where} Named set '{self.name}'")

                case 'X':
                    self.init_external(tokens[1:], where)

                case 'M':
                    for token in tokens[1:]:
                        self.independent_momenta.append(token)
                        self.independent_products += [(q,token) for q in self.independent_momenta]
                        logger.debug(f"{where} Added independent momentum '{token}'")

                case 'R':
                    if len(tokens) != 3:
                        raise ChPTError("Expected exactly two tokens (left- and right-hand side) in a replacement")
                    self.replacements.append( (tokens[1], Momentum(tokens[2], self.valid_momenta())) )
                    logger.debug(f"{where} Added replacement {tokens[1]} -> {tokens[2]}")

                case 'L':
                    for token in tokens[1:]:
                        self.loop_momenta.append(token)
                        logger.debug(f"{where} Added loop momentum '{token}'")

                case 'B':
                    expected_len = 3 + len(self.propagators + self.independent_products)
                    #if len(tokens) != expected_len:
                        #raise ChPTError(f"Invalid number of elements in propagator basis relation ({len(tokens)}, expected {expected_len})")

                    self.scalar_products[(tokens[1],tokens[2])] = tokens[3:]
                    logger.debug(f"""{where} Loaded relation {tokens[1]}*{tokens[2]} = {' '.join(
                            (f'- {coeff[1:]}*' if coeff.startswith('-') else f'+ {coeff}*' if coeff != '1' else '')
                            + (f'({self.propagators[i].momentum})^2'
                                if i < len(self.propagators) else
                                '.'.join(str(q) for q in self.independent_products[i - len(self.propagators)]))
                            for i, coeff in enumerate(tokens[3:]) if coeff != '0')
                        }""")

                case 'I':
                    for filename in tokens[1:]:
                        logger.debug(f"{where} Importing file '{filename}'...")
                        with open(filename, 'r') as file:
                            for n,l in enumerate(file):
                                diagram = self.parse_statement(f"{where}>{filename}:{n+1}", l, diagram)

                case 'P':
                    prop = Propagator(tokens[1:], self.valid_momenta())

                    if prop.momentum in self.propagator_name_map:
                        raise ChPTError(f"More than one propagator with momentum {prop.momentum}")
                    self.propagator_name_map[prop.momentum] = len(self.propagators)

                    self.propagators.append(prop)
                    logger.debug(f"{where} Added propagator {len(self.propagators)}: {prop}")

                case 'D':
                    self.add_diagram(diagram)
                    diagram = ChPTDiagram(tokens[1], sum(self.legs.values()))
                    if diagram.name in self.diagrams:
                        raise ChPTError(f"More than one diagram with name {diagram.name}")
                    logger.debug(f"{where} Initialized diagram '{diagram.name}'")

                case 'V':
                    if diagram is None:
                        raise ChPTError("Vertex must be inside diagram (expected 'D' before 'V')")

                    if len(tokens) >= 2 and tokens[1][0] == '@':
                        reference = tokens[1][1:]
                        if reference not in self.diagrams:
                            raise ChPTError(f"Unknown reference to diagram '{reference}'")

                        logger.debug(f"{where} Importing vertices from diagram '{reference}'...")
                        for vertex in self.diagrams[reference].vertices:
                            diagram.add_vertex(copy(vertex), where)
                    else:
                        diagram.add_vertex( Vertex(tokens[1:]), where )

                case 'E':
                    if diagram is None:
                        raise CHPTError("Edge must be inside diagram (expected 'D' before 'E')")

                    if len(tokens) >= 2 and tokens[1][0] == '@':
                        reference = tokens[1][1:]
                        if reference not in self.diagrams:
                            raise ChPTError(f"Unknown reference to diagram '{reference}'")

                        logger.debug(f"{where} Importing edges from diagram '{reference}'...")
                        for edge in self.diagrams[reference].edges:
                            diagram.add_edge(copy(edge), where)
                    else:
                        diagram.add_edge( Edge(diagram.vertices, diagram.vertex_name_map, self.propagators, self.propagator_name_map, tokens[1:]), where )

                case 'S':
                    if diagram is None:
                        raise ChPTError("Symmetry factor must be inside diagram (expected 'D' before 'S')")

                    if len(tokens) < 2:
                        raise ChPTError(f"Missing symmetry factor in S-statement")

                    diagram.set_symmetry_factor(tokens[1:])

                case 'G':
                    if diagram is None:
                        raise CHPTError("Group must be inside diagram (expected 'D' before 'G')")

                    if len(tokens) >= 2 and tokens[1][0] == '@':
                        reference = tokens[1][1:]
                        if reference not in self.diagrams:
                            raise ChPTError(f"Unknown reference to diagram '{reference}'")

                        logger.debug(f"{where} Importing symmetry group from diagram '{reference}'...")
                        diagram.permutation_group = self.diagrams[reference].permutation_group
                    else:
                        diagram.set_permutation_group(tokens[1:])

                case 'F':
                    if diagram is None:
                        raise ChPTError("Flags must be inside diagram (expected 'D' before 'F')")

                    diagram.add_flags(tokens[1:])

                case _:
                    raise ChPTError(f"Unknown statement '{tokens[0]}'")

        except Exception as err:
            errstr = str(err)
            if not errstr.endswith('\n'):
                errstr = errstr + '\n'
            raise ChPTError(errstr + indent(f"at {where}{newline(1)}{line}", ' '*4))

        return diagram

    def init_external(self, tokens, where=''):

        for token in tokens:
            if not specify_legs(token, self.legs):
                raise ChPTError("Invalid number-of-legs specification")

        logger.debug(f"{where} External legs: {'.'.join(f'{n}{p}' for p,n in self.legs.items())}")

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

        logger.debug(f"Finalized diagram '{diagram.name}'")

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
        print(dedent("""\
            #procedure trivialkinematics(DIAGRAM)
                .sort
                #call nskip(`DIAGRAM')"""),
            file=formfile)

        #Straightforward momentum relations
        for i,prop in enumerate(self.propagators):
            if prop.flav_dependent:
                continue
            print(indent(dedent(f"""\
                    id prop({+prop.momentum.substituted(self.replacements)},{prop.mass_squared}) = prop{i+1};
                    id prop({-prop.momentum.substituted(self.replacements)},{prop.mass_squared}) = prop{i+1};"""), " "*4),
                    #id [{prop.momentum}].[{prop.momentum}] = 1/prop{i+1} + {prop.mass_squared};
                    file=formfile)
        #for prop in self.all_propagators():
            #print(f"    id [{prop.momentum}] = {prop.momentum};", file=formfile)

        print(dedent(f"""\
                        .sort:>>trivial kinematics<<;
                    #endprocedure

                    #procedure fullkinematics(DIAGRAM)
                        #call trivialkinematics(`DIAGRAM')
                        #call nskip(`DIAGRAM')
                    """), file=formfile)

        # Basis relations
        for (p,q), coeffs in self.scalar_products.items():
            print(indent(dedent(f"""\
                    id {p}.{q} = {' '.join(
                        (f'- {coeff[1:]}*' if coeff.startswith('-') else f'+ {coeff}*' if coeff != '1' else '')
                        + (f'(1/prop{i+1} + {self.propagators[i].mass_squared})'
                           if i < len(self.propagators) else
                           f'({".".join(str(q) for q in self.independent_products[i - len(self.propagators)])})')
                        for i, coeff in enumerate(coeffs) if coeff != '0')
                    };"""), " "*4),
                    file=formfile)

        print(dedent(f"""\
                    .sort:>>full kinematics<<;
                #endprocedure

                #procedure replacements(DIAGRAM)
                    .sort
                    #call nskip(`DIAGRAM')

                    repeat;
                        #call replacementlist
                        #ifdef `REPLARG'
                            argument;
                                #call replacementlist
                            endargument;
                        #endif
                    endrepeat;
                #endprocedure

                #procedure replacementlist()
                    {newline(5).join(f'id {lhs} = {rhs};' for lhs,rhs in self.replacements)}
                #endprocedure
                """),
                file=formfile)

    def sp_symbols(self):
        return {s : s for s in self.external_momenta() + self.loop_momenta + self.independent_momenta}

    def print_header_FORM(self, formfile, options):
        print("#-", file=formfile)
        print(f"#appendpath {os.getcwd()}", file=formfile)

        for opt in options:
            if opt[0] == '-':
                continue

            eq = opt.find('=')
            if eq >= 0:
                if opt[:eq] == 'NF':
                    Nf = int(opt[eq+1:])
                print(f'#define {opt[:eq]} "{opt[eq+1:]}"', file=formfile)
            else:
                if opt == 'KALPAR':
                    Nf = "2"
                print(f"#define {opt}", file=formfile)

        print(dedent(f'''\
                #define NEXT "{sum(self.legs.values())}"
                #define ORDER "{order_OtoN(self.order)}"
                #define LOOPMOMENTA "{",".join(self.loop_momenta)}"
                #define INTMOMENTA "{",".join(self.independent_momenta)}"
                #include- ChPTdiagram.hf'''), file=formfile)

    def finalize_diagrams_FORM(self, vertex_roster, formfile):
        print(dedent(f"""\
                .sort:>>all {self.name} diagrams<<;
                off statistics;
                drop {', '.join(f'vert{vertex.name_FORM(tag)}' for vertex,count in vertex_roster.items() for tag in range(count))};
                """), file=formfile)

    def check_CoM(self):
        print("Checking conservation of momentum...")
        conservation = True
        for diagram in self.diagrams.values()   :
            violations = diagram.check_CoM(self.propagators, self.sp_symbols(), self.replacements)
            for vertex,perm,total in violations:
                print(f"    CoM violation in diagram {diagram.name}, vertex {vertex}, permutation {perm}: momenta sum to {total}")
                conservation = False

        if conservation:
            print("    CoM satisfied in all diagrams.")

        #filename = f"ChPTdiagram_{self.name}_CoM.frm"
        #with open(filename, 'w') as formfile:
            #print("off statistics;", file=formfile)
            #self.print_definitions_FORM(formfile)

            #for diagram in self.diagrams.values():
                #diagram.print_CoM(self.propagators, formfile)

            #self.print_replacements_FORM(formfile)

            #print(dedent("""\
                #.sort
                ##define RESULT "0"
                ##do EXPR={`ACTIVEEXPRNAMES_'}
                    ##if (termsin(`EXPR') > 0)
                        ##write " [FORM] diagram:vertex `EXPR' violates CoM by %E",`EXPR'
                        ##redefine RESULT "{`RESULT'+1}"
                    ##endif
                ##enddo
                ##if `RESULT' == 0
                    ##write " [FORM] Conservation of momentum is satisfied.\\n"
                ##else
                    ##write "\\n [FORM] Conservation of momentum failed in `RESULT' places.\\n"
                ##endif
                ##terminate `RESULT'
                #.end
                #"""), file=formfile)

        #print(dedent("""\
            #Running FORM to check conservation of momentum...
            #"""))
        #try:
            #subprocess.run(["form", "-q" ,filename], capture_output=False)
        #except subprocess.SubprocessError as err:
            #print(dedent(f'''\
                  #Running FORM failed; reason stated: "{err}"
                  #You can try running form -q {filename} manually instead.
                  #'''))

    def draw_graphs(self):
        filename = f"ChPTdiagram_{self.name}_graphs.m"

        with open(filename, 'w') as mathfile:
            print(dedent(f"""\
                (* ChPTdiagram, by Mattias Sjö 2024
                 *
                 * Mathematica file generated by ChPTdiagram.py on {datetime.now().strftime("%c")}
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

        print(dedent(f"""
              Mathematica code for drawing Feynman graphs has been written to {filename}.
              It is given as a (diagram name) -> (graph) association named '{self.name}graphs'
               and can be imported into a notebook for viewing.
              """))

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

        print(dedent(f"""\
            The following is the projection of scalar products onto the basis of inverse propagators.
            Following each 'B' are the two vectors in the product and then {len(vector_prods)} numbers.
            The ith number represents the coefficient in front of the squared momentum of propagator i
            (i.e., the inverse propagator plus the squared mass) in the scalar product.
            This can be pasted into the input file for future use."""))
        for elem in basis:
            print(elem)

    # No longer necessary: ChPTdiagram_vertex generates vertices on the fly
    #def ensure_vertices(self, vertex_roster, Nf, options):
        #for vertex in vertex_roster:
            #directory = Path(__file__).resolve().parent
            #vertex_path = Path(
                      #f"{directory}/vertices/"
                    #+ f"{''.join(f'{p}{vertex.legs.get(p,0)}' for p in PARTICLES)}"
                    #+ f"P{vertex.order}{f'NF{Nf}' if Nf else 'NF2' if 'KALPAR' in options else ''}"
                    #+ f"{''.join(f'_{opt}' for opt in options if opt[0] != '-')}.hf")
            #logger.debug(f"Vertex requested from {vertex_path}")
            #if not vertex_path.is_file():
                #logger.info(
                      #f"Generating missing O(p^{vertex.order}) "
                    #+ f"{' '.join(f'{n}-{PARTICLES[p]}' for p,n in vertex.legs.items())} "
                    #+ f"vertex {f'at {Nf=} ' if Nf else ''}"
                    #+ (f"with options {' '.join(options)}"  if options else ''))
                #FORM_command = (
                      #f"tform -l -q{''.join(f' -d N{p}={vertex.legs.get(p,0)}' for p in PARTICLES)}"
                    #+ f" -d NP={vertex.order}{f' -d NF={Nf} ' if Nf else ' '}"
                    #+ f"{' '.join(opt if opt[0] == '-' else f'-d {opt}' for opt in options)} "
                    #+  " ChPTdiagram_lagrangian.frm")
                #logger.debug(f"Vertex does not exist. Generating with '{FORM_command}'")
                #generate_vertex = subprocess.run(FORM_command.split(), cwd=str(directory))
                #if generate_vertex.returncode:
                    #raise ChPTError(f"Failed to generate missing vertex: FORM returned {generate_vertex.returncode}")
                #if not vertex_path.is_file():
                    #raise ChPTError(f"Failed to generate missing vertex: {vertex_path} not created")




    def generate_FORM(self, options = []):
        dirname = f"ChPTdiagram_{self.name}"
        os.makedirs(f"./{dirname}", exist_ok=True)
        os.makedirs(f"./{dirname}/diagrams", exist_ok=True)
        os.makedirs(f"./{dirname}/loops", exist_ok=True)
        os.makedirs(f"./{dirname}/flags", exist_ok=True)
        os.makedirs(f"./{dirname}/permute", exist_ok=True)

        Nf = None
        vertex_roster = {}

        # Have each diagram write its definitions, accumulating a total roster of vertices needed
        diagram_names = [
            diagram.define_FORM(vertex_roster, self.propagators, self.legs, set(self.loop_momenta), dirname, prefix = f"{self.name}x")
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

        with open(f"{dirname}/definitions.hf", 'w') as formfile:
            print_info_FORM(formfile, "This file defines all variables specific to this calculation.")

            self.print_definitions_FORM(formfile)

            print(f'#define DIAGRAMNAMES "{",".join(diagram_names)}"', file=formfile)
            print(f'#define VERTEXNAMES "{",".join(vertex_names)}"', file=formfile)

        with open(f"{dirname}/kinematics.hf", 'w') as formfile:
            print_info_FORM(formfile, "This file defines procedures for simplifying the kinematics.")
            self.print_kinematics_FORM(formfile)

    def generate_FORM_main(self, options = []):
        dirname = f"ChPTdiagram_{self.name}"
        filename = f"{dirname}.frm"

        self.generate_FORM(options)

        with open(filename, 'w') as formfile:

            # Creation message
            print_info_FORM(formfile, dedent(f"""\

                This is the main FORM file.
                Run with form -l {filename} (or tform -wNTHREAD to taste).
                Results are written to 'save/{self.name}.sav'.
                """), auto=False)

            # Define options, include relevant headers
            self.print_header_FORM(formfile, options)

            print(dedent(f'''\
                #include- {dirname}/definitions.hf
                #include- {dirname}/kinematics.hf

                * If diagrams are split into multiple expressions (say, projected onto different tensor structures),
                *  modify this so it nskips all variants of DIAGRAM's name.
                #procedure nskip(DIAGRAM)
                    skip; nskip `DIAGRAM';
                #endprocedure

                * This creates all the vertices, and may be a bit time-consuming.
                * The extraction from the Lagrangian is cached, so it should run faster in subsequent runs.
                #include- {dirname}/vertices.hf

                * Main loop - this will create and process each diagram in turn
                * Each repetition is done with the active diagram's flags enabled.
                #do DIAGRAM={{`DIAGRAMNAMES'}}
                    .sort
                    skip;
                    #include- {dirname}/flags/on_`DIAGRAM'.hf

                *    NOTE: any custom manipulations inserted between either of these steps should
                *    be guarded with #call nskip(`DIAGRAM')

                *    This sets up the Feynman rules of the diagram and contracts the flavor indices
                    #include- {dirname}/diagrams/`DIAGRAM'.hf

                *    This performs all permutations of the diagram's external legs
                    #include- {dirname}/permute/`DIAGRAM'.hf

                *    This enacts all replacements specified with 'R'
                *    It only applies inside function arguments if REPLARG is defined
                    #call replacements(`DIAGRAM')

                *    This identifies loop integrals
                    #call fullkinematics(`DIAGRAM')
                    #include- {dirname}/loops/`DIAGRAM'.hf{f'''
                    if(match(prop(?a)));
                        print "NOTE: in `DIAGRAM': %t";
                        #ifndef `REPLARG'
                            print "NOTE: REPLARG not defined, which might cause this";
                        #endif
                        exit "ERROR: failed to substitute all propagators";
                    endif;'''
                    if not any(p.flav_dependent for p in self.all_propagators()) else ''}

                    #include- {dirname}/flags/off_`DIAGRAM'.hf
                #enddo

                .sort:>>all {self.name} diagrams<<;
                off statistics;
                drop `VERTEXNAMES';

                * Room for custom post-processing.

                * Print and store results
                bracket `LECS', F, i_;
                print +s;
                .sort
                #call store({self.name})
                .end'''), file=formfile)

        # Ensure all necessary vertices are actually precomputed
        # No longer necessary: ChPTdiagram_vertex generates vertices on the fly
        #self.ensure_vertices(vertex_roster, Nf, options)

        #logger.info("User-defined #procedure projectexternal(DIAGRAM) needed to project out external flavors and polarizations")

        print(dedent(f"""
              FORM code for computing the diagrams has been written to {filename}.
              Run with form -l {filename} (or tform -wNTHREAD to taste).
              """))

def print_info_FORM(formfile, extra=[], auto=True):
    print(dedent(f"""\
        *** ChPTdiagram, by Mattias Sjö 2024
        *** FORM file generated by ChPTdiagram.py on {datetime.now().strftime("%c")}
        ***  with arguments '{' '.join(sys.argv[1:])}'"""), file=formfile)
    if auto:
        print(dedent(f"""\
            *** This file is automatically generated and should preferably not be modified.
            *** Instead, the main FORM file is meant to accommodate customization."""), file=formfile)
    for line in extra.split('\n'):
        print(f"*** {line}", file=formfile)

    print("", file=formfile)

def print_help():
    print(dedent('''\
        ChPTdiagram.py, by Mattias Sjö 2024

        The first argument is the name of an .xpt file containing the definition of a set of diagrams.
        Subsequent arguments are one or more of the following commands:
        - Help: print this help message and exit.
        - XPT: explain the syntax of the input .xpt file and exit.
        - CoM: check that conservation of momentum is respected by all input diagrams
        - Basis: determine how any scalar product of loop and external momenta can be expressed
            in terms of the basis of inverse propagators provided.
            This is both loaded into the program and printed as a set of B, which can be
            pasted into the .xpt file for future use.
        - Graphs: produce Mathematica code for drawing (very ugly) graphs of the diagrams,
            for verifying that the topologies have been properly implemented.
        - FORM: produce FORM code for computing the diagrams.
            Subsequent arguments are not parsed as commands, but are interpreted as FORM
            definitions, equivalent to what would be passed into FORM itself using -d flags.
        - FORM-main: generate FORM files as above, and produce a template FORM main file set up to include
            those files and perform the actual calculation, while providing convenient spaces to customize
            the computations.
            This is a separate command from FORM to avoid overwriting such customizations when regenerating
            the FORM files.
        Commands are case-insensitive, and may be abbreviated.
        '''))
def print_xpt_help():
    print(dedent('''\
        ChPT file format, by Mattias Sjö 2024

        This is a very simple way of describing multiloop diagrams in ChPT
         using a formulation suitable for master integral reduction.
        All non-empty lines contain a single statement, consisting of a single character
         specifying the type of statement followed by a number of whitespace-separated
         arguments.
        The statement types are the following:
        - #: comment, ignore the rest of the line.
        - I: input the contents of other .xpt files given as arguments.
        - N: name the diagram set according to the single argument.
             Exactly one N-statement must be given.
        - X: Specify the number of external legs (all diagrams have the same).
             Each argument is a number followed by a character specifying the type of
             particle. The particle types are M(eson), V(ector), A(xial vector), S(calar) and
             P(seudoscalar). The particle specifier defaults to M(eson) if omitted. The
             number of each type of particle may be specified only once, and defaults to
             zero. External legs are assigned incoming momenta pi, with i being a 1-based
             index, assigned incrementally with mesons first, then vectors, etc.
        - L: Define loop momenta, each one given as a separate argument.
             Bear in mind that l1, l2, etc. may clash with the names of the LECs.
        - M: Define external momenta, each one given as a separate argument.
             All independent external momenta should be explicitly given this way.
        - R: Specify a replacement to be applied to the expressions.
             For example, this may be used to implement conservation of external momentum, or
             to replace the automatically assigned external momenta with user-defined ones.
             Takes two arguments, the left- and right-hand side of the replacement,
             respectively. Example: "R p4 -(p1+p2+p3)" for 4-point diagrams. Bear in mind
             that the expressions should work in both FORM and Mathematica.
        - P: Define a propagator.
             Takes two arguments: the momentum and mass-squared of the propagator.
             Propagators are numbered in the order they are specified, starting at 1. Bear in
             mind that for master integral reduction, all scalar products of loop momenta,
             and all scalar products of a loop momentum and an independent external momentum,
             must be expressible as a linear combination of squared propagator momenta. Thus,
             the number of propagators should be E*L + L*(L+1)/2 with L loops and E
             independent external momenta; this may require the inclusion of propagators that
             don't actually appear as a propagator in any diagram. The command 'basis' checks
             that this has been done correctly. The mass-squared may be given as '*',
             indicating flavor-dependent mass.
        - B: Specify how a product of momenta is expressed as a linear combination of squared
             propagator momenta, as described under 'P'.
             The first two arguments are the momenta being multiplied together, and the
             following are the coefficients in front of each squared propagator momentum, in
             the order the propagators were defined. The command 'basis' automatically
             generates all necessary B-statements.
        - D: Initialize a new active diagram, and finalize the previous one, if any.
             [VESGF]-statements act on an active diagram, so a D-statement must be made
             before these can be used. It is recommended to give all other types of
             statements before the first D-statement, as a preamble. The first argument is
             compulsory, and gives the name of the diagram; this is combined with the diagram
             set name to produce the full diagram name. The second, optional argument gives
             the symmetry factor of the diagram. Exponents (^) and factorials (!) may be used
             and will automatically be translated to FORM and Mathematica syntax.
        - V: Define a vertex for the active diagram.
             There are three types of arguments, and they can be given in any order. If the
             argument is of the form of zero or more N's followed by 'LO', or alternatively
             'N', a number representing the number of N's, and 'LO', that specifies the order
             of the vertex, defaulting to LO if omitted. If the argument would be accepted by
             an X-statement, it specifies the number of legs on the vertex as described
             there. Otherwise, the argument specifies the optional name of the vertex.
             Alternatively, there may be a single argument, consisting of '@' and the name of
             a previously defined diagram; this imports all vertices of that diagram as if
             all its V-statements were pasted in place of the 'V @...' Vertices are numbered
             in the order they are defined, starting at 1.
        - E: Define an edge for the active diagram.
             There are three arguments: the number or name of the source vertex, the number
             or name of the destination vertex, and the number or momentum of the propagator
             represented by the edge. The momentum of the propagator flows from the source to
              the destination. The propagator may instead be given as '*MASS', in which case
             the momentum will be deduced from conservation. MASS should be the mass-squared
             of the propagator, or another '*' to indicate flavor-dependent mass.
             Alternatively, there may be a single argument, consisting of '@' and the name of
             a previously defined diagram; this imports all edges of that diagram similarly
             to 'V @...'. Note that vertex names and numbers may be different in this
             diagram, so this should be used with caution.
        - S: Endow the current diagram with a symmetry factor.
             The first argument should be an expression for the symmetry factor, by which the
             diagram will be divided. Exponents (^) and factorials (!) may be used and will
             automatically be translated to FORM and Mathematica syntax.
             Optional symmetry factors may be given, prefixed by the name of a FORM
             preprocessor variable and a colon. If that variable is defined at the moment the
             diagram is created, that factor will be used instead (the first one in the order
             given here takes precedence). There is an implicit NOSYMFACT:1 that overrides
             other alternative symmetry factors.
        - G: Endow the current diagram with a group of permutations of its external legs,
             over which it should be summed. There can be at most one G-statement per diagram.
             The default arrangement of external momenta is that all momenta belonging to
             each particle type are assigned in the order they appeared in the X-statements,
             and for each type, they are assigned to all available legs on each vertex in the
             order they appeared in the V-statements. The group can be specified in either of
             three ways: one or more permutations, a distribution, or as a copy of another
             diagram's group with 'G @...'. Permutations are given in cycle notation
             [example: (1,2)(5,6,7,8)], and the group will be all distinct permutations that
             can be formed by composing the specified permutations (including the identity).
             A distribution is given as a bracket-enclosed list (example: [2,2,4]) whose
             elements sum to the number of legs of the diagram. The group will be all ways of
             distributing the external momenta into sets with the given list of sizes,
             without regard to the order within each set, and without regard to the relative
             order of sets of equal size. Each set-size may be optionally suffixed by an
             alphanumeric label, and in that case the relative order is respected between
             sets with different labels; the special label '*' is considered different from
             all labels including itself. Normally, each set will correspond to the number of
             external legs on each vertex, and the group will be all ways of assigning the
             external legs, up to symmetries of the diagram which can be accounted for with
             the labels; for instance, [2,2] is appropriate for a simple lollipop diagram
             with four identical external legs, two on each vertex; [2a,2b] or [2*,2*] is
             appropriate if the vertices are not interchangeable.
        - F: Equip a diagram with one or more flags.
             By "flags" is meant FORM preprocessor variables that will be defined only while
             processing the current diagram, for convenience when controlling any diagram-
             specific procedures. The syntax is as for FORM -d options, i.e. NAME or NAME=VALUE.
'''))

def main():

    logging.basicConfig(level=logging.DEBUG)

    if len(sys.argv) < 2:
        print("ERROR: No arguments given. Write argument 'help' for a list of valid arguments.")
        sys.exit(1)

    if 'help'.startswith(sys.argv[1].casefold()):
        print_help()
        return
    if 'xpt'.startswith(sys.argv[1].casefold()):
        print_xpt_help()
        return

    try:
        diagrs = ChPTDiagramSet(sys.argv[1])
    except FileNotFoundError as err:
        print(f"ERROR: xpt file not found: '{sys.argv[1]}'", file=sys.stderr)
        sys.exit(1)
    except ChPTError as err:
        print(f"ERROR: {err}", file=sys.stderr)
        sys.exit(1)

    for i,arg in enumerate(sys.argv[2:]):
        if 'help'.startswith(arg.casefold()):
            print_help()
            break
        if 'xpt'.startswith(arg.casefold()):
            print_xpt_help()
            break
        elif 'com'.startswith(arg.casefold()):
            diagrs.check_CoM()
            continue
        elif 'graphs'.startswith(arg.casefold()):
            diagrs.draw_graphs()
            continue
        elif 'basis'.startswith(arg.casefold()):
            diagrs.determine_prop_basis()
            continue
        elif 'form'.startswith(arg.casefold()):
            diagrs.generate_FORM(options=sys.argv[2:][i+1:])
            break
        elif 'form-main'.startswith(arg.casefold()):
            diagrs.generate_FORM(options=sys.argv[2:][i+1:])
            diagrs.generate_FORM_main(options=sys.argv[2:][i+1:])
            break
        #elif 'symmetry'.startswith(arg.casefold()):
            #diagrs.determine_symmetry(options=sys.argv[2:][i+1:])
            #break
        else:
            raise ChPTError(f"Unrecognized command: {arg} (ChPTdiagram help' gives a list of valid ones)")

if __name__ == '__main__':
    main()
