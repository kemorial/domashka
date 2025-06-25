import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from calc_pack import run_simulation, Data

def plot_results(results):
    fig = plt.figure(figsize=(25, 25))
    gs = GridSpec(ncols=3, nrows=3, figure=fig, wspace=0.5, hspace=1)
    
    # Графики параметров
    plots = [
        (gs[0, 0], 'phi', results.phi, 'Угол поворота'),
        (gs[0, 1], 'dL', results.dL, 'Элементарная работа'),
        (gs[0, 2], 'dV', results.V, 'Объём'),
        (gs[1, 0], 'p', results.p, 'Давление', results.p1, 'Без сгорания'),
        (gs[1, 1], 'dQx', results.dQx, 'Тепловыделение'),
        (gs[1, 2], 'cv', results.cv[:-1], 'Теплоёмкость'),
        (gs[2, 0], 'T', results.T[:-1], 'Температура'),
        (gs[2, 1], 'dQw', results.dQw[:-1], 'Теплопотери'),
        (gs[2, 2], 'pV', (results.V, results.p), 'Индикаторная диаграмма')
    ]
    
    for pos, ylabel, data, title, *compare in plots:
        ax = plt.subplot(pos)
        if ylabel == 'pV':
            ax.plot(*data)
            ax.set_xlabel("V")
            ax.set_ylabel("p")
        else:
            ax.plot(results.shk[:len(data)], data)
            ax.set_ylabel(ylabel)
            if compare:
                ax.plot(results.shk[:len(compare[0])], compare[0])
                ax.legend([title, compare[1]])
        ax.set_title(title)
        ax.grid()
    
    # Дополнительные графики
    fig2 = plt.figure(figsize=(25, 25))
    gs2 = GridSpec(ncols=1, nrows=3, figure=fig2, wspace=1, hspace=1)
    
    plots2 = [
        (gs2[0], results.phi, results.p, 'Давление от угла поворота', 'phi', 'P'),
        (gs2[1], results.phi, results.T, 'Температура от угла поворота', 'phi', 'T'),
        (gs2[2], results.V, results.p, 'Индикаторная диаграмма', 'V', 'P')
    ]
    
    for pos, x, y, title, xlabel, ylabel in plots2:
        ax = plt.subplot(pos)
        ax.plot(x, y)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid()
    
    # Вывод параметров
    print(f"Индикаторное давление: {results.p_ind:.2f} Па")
    print(f"Индикаторная мощность: {results.N_ind:.2f} Вт")
    print(f"Удельный индикаторный расход: {results.g_ind:.5f} кг/Дж")
    print(f"Индикаторный КПД: {results.ef_ind:.3f}")
    print(f"Работа цикла: {results.L:.2f} Дж")
    
    plt.show()

if __name__ == "__main__":
    data = Data()
    results = run_simulation(data)
    plot_results(results)