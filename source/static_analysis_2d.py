import os
import sympy
import networkx
import numpy
import json
import matplotlib.pyplot as plt


def main(
    liaisons: dict[str, list[int | float, int | float, bool, bool, bool]],
    parts: dict[str, list[str]],
    forces: dict[str, list[int | float, int | float, int | float, int | float, int | float]],
) -> dict[str, int | float]:
    """
    Analyse statique 2D des efforts dans les liaisons d'une structure

    :param liaisons: link_name[x, y, transmit_fx, transmit_fy, transmit_rz]
    :param parts: part_name[link_name_1, link_name_2, ...]
    :param forces: part_name[[pos_f_x_1, pos_f_y_1, fx_1, fy_1, rz_1], [pos_f_x_2, pos_f_y_2, fx_2, fy_2, rz_2], ...]
    :return:
    """
    check_inputs(liaisons, parts, forces)
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
    solution = solve(equations, symbols)
    result = {str(key): val for key, val in solution.items()} if solution else {}
    print("\nRESULTS:\n")
    for force, value in result.items():
        if value == 0:
            continue
        print(f"{force} = {round(value, 10)}")
    # display_liaison_graph(liaisons, parts, forces, result, liaison_parts)
    return result


def check_inputs(
    liaisons: dict[str, list[int | float, int | float, bool, bool, bool]],
    parts: dict[str, list[str]],
    forces: dict[str, list[int | float, int | float, int | float, int | float, int | float]],
) -> None:
    if "GND" not in parts:
        raise Exception(f"(parts): Missing 'GND' part in {parts}")
    for liaison in liaisons:
        counter = 0
        for part, liaions_part in parts.items():
            for liaion_part in liaions_part:
                if liaion_part not in liaisons:
                    raise Exception(f"(liaisons): Link {liaion_part} not referenced in {liaisons}")
            if liaison in liaions_part:
                counter += 1
        if counter != 2:
            raise Exception(f"(liaisons): Liaison {liaison} only into {counter} part, has to be 2")
    for part_force in forces:
        if part_force not in parts:
            raise Exception(f"(forces): Part {part_force} missing in {parts}")


def hyperstatisme(
    symbols_indexes: dict[str, int],
    symbols_null: list[int],
    parts: dict[str, list[str]],
) -> int:
    ns = len(symbols_indexes) - len(symbols_null)
    p = len(parts)
    hyperstatisme_degree = ns - 3 * (p - 1)
    print(f"Hyperstatisme = {hyperstatisme_degree}")
    if hyperstatisme_degree != 0:
        raise Exception(f"Mecanism not Isostatic (degree {hyperstatisme_degree})\n\t{ns = }\n\t{p = }")
    return hyperstatisme_degree


def get_symbols(
    liaisons: dict[str, list[int | float, int | float, bool, bool, bool]],
) -> tuple[sympy.symbols, dict[str, int]]:
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
    return (
        sympy.symbols(symbols),
        symbols_indexes
    )


def get_symbols_null(
    liaisons: dict[str, list[int | float, int | float, bool, bool, bool]],
) -> list[int]:
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


def get_equations_null(
    symbols: sympy.symbols,
    symbols_null: list[int],
) -> list[sympy.Eq]:
    equations = []
    for symbol_null in symbols_null:
        equation = sympy.Eq(symbols[symbol_null], 0)
        equations.append(equation)
    return equations


def get_equations_fx(
    symbols: sympy.symbols,
    parts: dict[str, list[str]],
    symbols_indexes: dict[str, int],
    forces: dict[str, list[int | float, int | float, int | float, int | float, int | float]],
    symbols_null: list[int],
) -> list[sympy.Eq]:
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
        equation = sympy.Eq(sum(symbols_sum), 0)
        print(f"\t{part} - {equation}")
        equations.append(equation)
    return equations


def get_equations_fy(
    symbols: sympy.symbols,
    parts: dict[str, list[str]],
    symbols_indexes: dict[str, int],
    forces: dict[str, list[int | float, int | float, int | float, int | float, int | float]],
    symbols_null: list[int],
) -> list[sympy.Eq]:
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
        equation = sympy.Eq(sum(symbols_sum), 0)
        print(f"\t{part} - {equation}")
        equations.append(equation)
    return equations


def get_equations_rz(
    symbols: sympy.symbols,
    parts: dict[str, list[str]],
    symbols_indexes: dict[str, int],
    forces: dict[str, list[int | float, int | float, int | float, int | float, int | float]],
    symbols_null: list[int],
    liaisons: dict[str, list[int | float, int | float, bool, bool, bool]],
) -> list[sympy.Eq]:
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
            equation = sympy.Eq(sum(symbols_sum), 0)
            print(f"\t{part} - {liaison} - {equation}")
            equations.append(equation)
    return equations


def get_torque(
    fx: int | float,
    fy: int | float,
    rz: int | float,
    pfx: int | float,
    pfy: int | float,
    x: int | float,
    y: int | float,
) -> float:
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


def solve(equations, symbols):
    solution = sympy.nsolve(equations, symbols, [0 for _ in symbols], dict=True)[0]
    return solution


if __name__ == "__main__":
    main(**json.load(open(os.path.join(
        os.path.dirname(__file__), "..", "data",
        "chariot_droite.json"
    ))))
