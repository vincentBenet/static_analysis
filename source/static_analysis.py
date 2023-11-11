from math import atan, cos, sin, pi
import sympy as sym
import networkx
import numpy
import json
import matplotlib.pyplot as plt


def main(liaisons, parts, forces):
    check(liaisons, parts, forces)
    liaison_parts = get_liaison_parts(liaisons, parts)
    symbols, symbols_indexes = get_symbols(liaisons)
    symbols_null = get_symbols_null(liaisons)
    print(f"NB inconnues = {len(liaisons) * 3 - len(symbols_null)}")
    hyperstatisme(symbols_indexes, symbols_null, parts)
    equations_null = get_equations_null(symbols, symbols_null)
    equations_fx = get_equations_fx(symbols, parts, symbols_indexes, forces, symbols_null)
    equations_fy = get_equations_fy(symbols, parts, symbols_indexes, forces, symbols_null)
    equations_rz = get_equations_rz(symbols, parts, symbols_indexes, forces, symbols_null, liaisons)
    equations = equations_null + equations_fx + equations_fy + equations_rz
    print(f"NB equations = {len(equations) - len(equations_null)}")
    solution = sym.solve(equations, symbols, warn=True)
    result = {str(key): val for key, val in solution.items()} if solution else {}
    print("\nRESULTS:\n")
    for force, value in result.items():
        if value == 0:
            continue
        print(f"{force} = {value}")
    display_liaison_graph(liaisons, parts, forces, result, liaison_parts)
    return result


def get_liaison_parts(liaisons, parts):
    res = {}
    for liaison in liaisons:
        res[liaison] = []
        for part, liaisions_part in parts.items():
            if liaison in liaisions_part:
                res[liaison].append(part)
    return res

def check(liaisons, parts, forces):
    check_liaisons_types(liaisons)
    check_parts_types(parts)
    check_forces_types(forces)
    assert "GND" in parts
    for liaison in liaisons:
        counter = 0
        for part, liaions_part in parts.items():
            for liaion_part in liaions_part:
                assert liaion_part in liaisons
            if liaison in liaions_part:
                counter += 1
        if counter != 2:
            raise Exception(f"Liaison {liaison} only into {counter} part")
    for part_force in forces:
        assert part_force in parts


def check_liaisons_types(liaisons):
    for liaison, elems in liaisons.items():
        x, y, fx, fy, rz = elems
        assert isinstance(liaison, str)
        assert isinstance(x, int | float)
        assert isinstance(y, int | float)
        assert isinstance(fx, int | float)
        assert isinstance(fy, int | float)
        assert isinstance(rz, int | float)


def check_parts_types(parts):
    for part, liaisons in parts.items():
        assert isinstance(part, str)
        for liaison in liaisons:
            assert isinstance(liaison, str)


def check_forces_types(forces):
    for part, forces in forces.items():
        assert isinstance(part, str)
        for force in forces:
            x, y, fx, fy, rz = force
            assert isinstance(x, int | float)
            assert isinstance(y, int | float)
            assert isinstance(fx, int | float)
            assert isinstance(fy, int | float)
            assert isinstance(rz, int | float)


def hyperstatisme(symbols_indexes, symbols_null, parts):
    ns = len(symbols_indexes) - len(symbols_null)
    p = len(parts)
    hyperstatisme_degree = ns - 3 * (p - 1)
    print(f"Hyperstatisme = {hyperstatisme_degree}")
    if hyperstatisme_degree != 0:
        raise Exception(f"Mecanism not Isostatic (degree {hyperstatisme_degree})\n\t{ns = }\n\t{p = }")
    return hyperstatisme_degree


def get_symbols(liaisons):
    symbols = []
    symbols_indexes = {}
    i = 0
    for liaison in liaisons:
        fx = f"fx{liaison}"
        fy = f"fy{liaison}"
        rz = f"rz{liaison}"
        symbols_indexes[fx] = i
        symbols.append(fx)
        symbols_indexes[fy] = i + 1
        symbols.append(fy)
        symbols_indexes[rz] = i + 2
        symbols.append(rz)
        i += 3
    return sym.symbols(symbols), symbols_indexes


def get_symbols_null(liaisons):
    symbols_null = []
    i = 0
    for liaison, elems in liaisons.items():
        x, y, fx, fy, rz = elems
        if not fx:
            symbols_null.append(i)
        if not fy:
            symbols_null.append(i+1)
        if not rz:
            symbols_null.append(i+2)
        i += 3
    return symbols_null


def get_equations_null(symbols, symbols_null):
    equations = []
    for symbol_null in symbols_null:
        equation = sym.Eq(symbols[symbol_null], 0)
        equations.append(equation)
    return equations


def get_equations_fx(symbols, parts, symbols_indexes, forces, symbols_null):
    print("Equations Fx")
    equations = []
    for part, links in parts.items():
        if part == "GND":
            continue
        symbols_sum = []
        for liaison in links:
            symbol_str = f"fx{liaison}"
            index = symbols_indexes[symbol_str]
            if index not in symbols_null:
                symbol = symbols[index]
                symbols_sum.append(symbol)
        for force in forces.get(part, []):
            x, y, fx, fy, rz = force
            symbols_sum.append(fx)
        equation = sym.Eq(sum(symbols_sum), 0)
        print(f"\t{part} - {equation}")
        equations.append(equation)
    return equations


def get_equations_fy(symbols, parts, symbols_indexes, forces, symbols_null):
    print("Equations Fy")
    equations = []
    for part, links in parts.items():
        if part == "GND":
            continue
        symbols_sum = []
        for liaison in links:
            symbol_str = f"fy{liaison}"
            index = symbols_indexes[symbol_str]
            if index not in symbols_null:
                symbol = symbols[index]
                symbols_sum.append(symbol)
        for force in forces.get(part, []):
            x, y, fx, fy, rz = force
            symbols_sum.append(fy)
        equation = sym.Eq(sum(symbols_sum), 0)
        print(f"\t{part} - {equation}")
        equations.append(equation)
    return equations


def get_equations_rz(symbols, parts, symbols_indexes, forces, symbols_null, liaisons):
    print("Equations Rz")
    equations = []
    for liaison, elems in liaisons.items():
        x, y, _, _, rz = elems
        for part, liaisons_parts in parts.items():
            if part == "GND":
                continue
            symbols_sum = []
            for force in forces.get(part, []):
                cx, cy, fx, fy, rz = force
                symbols_sum.append(get_torque(fx, fy, rz, cx, cy, x, y))
            if liaison not in liaisons_parts:
                continue
            for liaison_part in liaisons_parts:
                cx, cy, _, _, _ = liaisons[liaison_part]
                index_fx = symbols_indexes[f"fx{liaison_part}"]
                index_fy = symbols_indexes[f"fy{liaison_part}"]
                index_rz = symbols_indexes[f"rz{liaison_part}"]
                fx = symbols[index_fx] if index_fx not in symbols_null else 0
                fy = symbols[index_fy] if index_fy not in symbols_null else 0
                rz = symbols[index_rz] if index_rz not in symbols_null else 0
                symbols_sum.append(get_torque(fx, fy, rz, cx, cy, x, y))
            equation = sym.Eq(sum(symbols_sum), 0)
            print(f"\t{part} - {liaison} - {equation}")
            equations.append(equation)
    return equations


def get_torque(fx, fy, rz, pfx, pfy, x, y):
    dist_x = pfx - x
    dist_y = pfy - y
    rz_fx = - dist_y * fx
    rz_fy = dist_x * fy
    torque = 0
    if rz_fx:
        torque += rz_fx
    if rz_fy:
        torque += rz_fy
    if rz:
        torque += rz
    return torque


def display_liaison_graph(liaisons_dict, parts, forces, result, liaison_parts):
    graph = networkx.DiGraph()
    colors = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]]
    color_index = 0
    force_index = 0
    pos = {}
    for part, liaisons in parts.items():
        n = len(liaisons)
        for i in range(n):
            liaison_1 = liaisons[i]
            pos_x_1, pos_y_1, inc_fx_1, inc_fy_1, inc_rz_1 = liaisons_dict[liaison_1]
            liaison_name_1 = get_liaison_name(inc_fx_1, inc_fy_1, inc_rz_1)
            fx_1 = result.get('fx' + liaison_1, 0)
            fy_1 = result.get('fy' + liaison_1, 0)
            rz_1 = result.get('rz' + liaison_1, 0)
            txt_1 = get_node_txt(liaison_1, liaison_name_1, fx_1, fy_1, rz_1, liaison_parts[liaison_1])
            if liaison_1 not in graph.nodes():
                pos[liaison_1] = numpy.array([pos_x_1, pos_y_1])
                graph.add_node(
                    liaison_1,
                    type_liaison=liaison_name_1,
                    fx=fx_1,
                    fy=fy_1,
                    rz=rz_1,
                    x=pos_x_1,
                    y=pos_y_1,
                    label=txt_1
                )
        for i in range(n):
            liaison_1 = liaisons[i]
            liaison_2 = liaisons[(i + 1) % n]
            txt = f"{part}"
            color = colors[color_index % len(colors)]
            graph.add_edge(liaison_1, liaison_2, label=txt, color=color)
        color_index += 1
    networkx.draw(
        graph, pos,
        # edge_color=networkx.get_edge_attributes(graph, 'color'),
        labels=networkx.get_node_attributes(graph, 'label'),
        node_size=[10000 for _ in range(len(graph.nodes()))],
    )
    networkx.draw_networkx_edge_labels(graph, pos, edge_labels=networkx.get_edge_attributes(graph, 'label'))
    for part, forces_added in forces.items():
        for force in forces_added:
            force_index += 1
            pos_x, pos_y, fx, fy, rz = force
            txt = f"F{force_index} -> {part}"
            if fx:
                txt += f"Fx={round(fx, 2)}N"
            if fy:
                txt += f"\nFy={round(fy, 2)}N"
            if rz:
                txt += f"\nRz={round(rz, 2)}N.m"
            plt.arrow(pos_x, pos_y, fx/1000, fy/1000, width=(fx**2 + fy**2)**0.5/10000)
            plt.text(pos_x, pos_y, txt)
    plt.show()


def get_node_txt(liaison, liaison_name, fx, fy, rz, liaison_parts):
    txt = f"{liaison} {'-'.join(liaison_parts)}\n{liaison_name}"
    if fx:
        txt += f"\nFx={int(fx*10)/10} N"
    if fy:
        txt += f"\nFy={int(fy*10)/10} N"
    if rz:
        txt += f"\nRz={int(rz*10)/10} N.m"
    return txt


def get_liaison_name(inc_fx, inc_fy, inc_rz):
    if inc_fx and inc_fy and inc_rz:
        return "Liaison encastrement"
    elif inc_fx and not inc_fy and inc_rz:
        return "Liaison glissière Y"
    elif inc_fx and inc_fy and not inc_rz:
        return "Liaison pivot Z"
    elif not inc_fx and inc_fy and inc_rz:
        return "Liaison glissière X"
    elif not inc_fx and not inc_fy and inc_rz:
        return "Liaison plane"
    elif not inc_fx and inc_fy and not inc_rz:
        return "Liaison sphère cylindre X"
    elif inc_fx and not inc_fy and not inc_rz:
        return "Liaison sphère cylindre Y"
    else:
        raise Exception("Pas le liaison")


if __name__ == "__main__":
    main(**json.load(open("chariot_droite.json")))
