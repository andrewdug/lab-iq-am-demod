#!/usr/bin/env python3
"""Цифровой обработчик к занятию по фильтру Гильберта (п. 3-8), бригада M = 1.

Вход: результаты моделирования Qucs-S (qucsator)
  LR2k5_N1_hilbert_AC.dat - АЧХ, ФЧХ фильтров и элементов задержки;
  LR2k5_N1_hilbert_TR.dat - x0 = s(t); xd, xd14, xd18 = s'(t); y, yw, y14, y18 = ŝ(t).
Выход (--out): рисунки fig_*.png, results.json с численными результатами,
hilbert_iq.csv с отсчётами s, s', ŝ, I, Q, U итогового фильтра.

Запуск из корня репозитория:  python scripts/hilbert_iq.py --dat qucs --out report/figures
Нужны numpy и matplotlib.
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Исходные данные
M = 1                      # номер бригады
F0 = (400 + 30 * M) * 1e6  # несущая, Гц
FM = 3 * M * 1e6           # частота модуляции, Гц
U0 = 1.0                   # амплитуда несущей, В
MOD = 1.0                  # коэффициент модуляции
FS = 4 * F0                # частота дискретизации
DT = 1 / FS
T_AN = (1e-6, 2e-6)        # окно анализа: 3 периода Fм = 1720 отсчётов (шаг ДПФ 1 МГц)
T_SKIP = 0.05e-6           # заполнение линии задержки, с
HOD_FM_PERIODS = 1         # годограф: один период модуляции

VARIANTS = [
    dict(key="o10", name="порядок 10, без окна", ref="xd", out="y", ac="y", acref="xd", order=10, color="#d62728"),
    dict(key="w10", name="порядок 10, Хэмминг", ref="xd", out="yw", ac="yw", acref="xd", order=10, color="#1f77b4"),
    dict(key="w14", name="порядок 14, Хэмминг", ref="xd14", out="y14", ac="y14", acref="xd14", order=14, color="#ff6600"),
    dict(key="w18", name="порядок 18, Хэмминг", ref="xd18", out="y18", ac="y18", acref="xd18", order=18, color="#2ca02c"),
]
MAIN = "w14"               # итоговый фильтр


def coeffs(order, hamming):
    """Коэффициенты h(n), n = 1..N, из формулы презентации ex5 (нулевые при чётном n - n0)."""
    n_coef = order + 1
    n0 = order // 2 + 1
    n = np.arange(1, n_coef + 1)
    k = n - n0
    h = np.zeros(n_coef)
    odd = k % 2 != 0
    h[odd] = 2 / (np.pi * k[odd])
    if hamming:
        h = h * (0.54 - 0.46 * np.cos(2 * np.pi * (n - 1) / (n_coef - 1)))
    return h


def freq_resp(h, f):
    d = np.arange(len(h)) * DT
    return np.exp(-2j * np.pi * np.outer(np.atleast_1d(f), d)) @ h


_CPLX = re.compile(r"^([+-]?[\d.]+(?:e[+-]?\d+)?)([+-])j([\d.]+(?:e[+-]?\d+)?)$", re.I)


def read_dat(path, names):
    """Читает из набора данных Qucs только нужные переменные."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    out = {}
    for m in re.finditer(r"<(indep|dep) (\S+)[^>]*>\s*(.*?)\s*</\1>", text, re.S):
        name = m.group(2)
        if name not in names:
            continue
        tok = m.group(3).split()
        if "j" in tok[0]:
            parts = [_CPLX.match(x) for x in tok]
            out[name] = np.array([complex(float(p[1]), float(p[2] + p[3])) for p in parts])
        else:
            out[name] = np.array(tok, dtype=float)
    miss = set(names) - set(out)
    if miss:
        raise SystemExit(f"В {path} нет переменных: {sorted(miss)}")
    return out


def spectrum(x):
    return np.fft.fftshift(np.abs(np.fft.fft(x))) / len(x)


def db(x):
    return 20 * np.log10(np.maximum(x, 1e-12))


def rms(x):
    return np.sqrt(np.mean(np.square(x)))


def hilbert_ref(t, s):
    """Эталон H{s}(t - τ) через БПФ: на 2 мкс укладывается целое число периодов f0 и Fм."""
    n = len(t) - 1
    spec = np.fft.fft(s[:n])
    freq = np.fft.fftfreq(n, (t[-1] - t[0]) / n)
    return spec * (-1j * np.sign(freq)), freq


def measure_delay(y, spec, freq, steady, tau_max):
    """МНК-совмещение y(t) с K·H{s}(t - τ): сетка, затем золотое сечение."""
    y = y[steady]

    def fit(tau):
        h = np.fft.ifft(spec * np.exp(-2j * np.pi * freq * tau)).real
        h = np.append(h, h[0])[steady]
        k = (h @ y) / (h @ h)
        return np.sum((y - k * h) ** 2), k

    grid = np.linspace(0, tau_max, 241)
    i = int(np.argmin([fit(x)[0] for x in grid]))
    a, b = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    g = (np.sqrt(5) - 1) / 2
    for _ in range(60):
        c, d = b - g * (b - a), a + g * (b - a)
        if fit(c)[0] < fit(d)[0]:
            b = d
        else:
            a = c
    tau = (a + b) / 2
    return tau, fit(tau)[1]


def phase_at(a, b, tn, f):
    e = np.exp(-2j * np.pi * f * tn)
    return np.angle(np.sum(a * e) / np.sum(b * e), deg=True)


def style_axes(ax, xlabel, ylabel, title=None):
    ax.grid(True, color="#cccccc", linewidth=0.6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=11)


def new_fig(rows=1, h=3.0):
    fig, ax = plt.subplots(rows, 1, figsize=(8.2, h * rows), constrained_layout=True)
    return fig, ax


def plot_ac(ac, out):
    f = ac["acfrequency"] / 1e6
    fig, ax = new_fig(2, 3.1)
    ax[0].plot(f, np.abs(ac["xd.v"]), "k--", lw=1.2, label="элемент задержки 5Δt")
    for v in VARIANTS:
        ax[0].plot(f, np.abs(ac[v["ac"] + ".v"]), color=v["color"], lw=1.2, label=v["name"])
    ax[0].axvline(F0 / 1e6, color="#888888", ls=":", lw=1)
    style_axes(ax[0], "f, МГц", "|K|", "АЧХ в диапазоне 1 кГц … fд/2 = 860 МГц (пунктир: f0 = 430 МГц)")
    ax[0].set_xlim(0, 860)
    ax[0].legend(fontsize=8, loc="lower center", ncol=2)
    sel = (f > 380) & (f < 480)
    for v in VARIANTS:
        ax[1].plot(f[sel], np.abs(ac[v["ac"] + ".v"])[sel], color=v["color"], lw=1.4, label=v["name"])
    ax[1].axhline(1.01, color="#888888", ls="--", lw=0.9)
    ax[1].axhline(0.99, color="#888888", ls="--", lw=0.9, label="допуск ±1 %")
    ax[1].axvspan(427, 433, color="#dddddd", label="полоса сигнала 427…433 МГц")
    style_axes(ax[1], "f, МГц", "|K|", "АЧХ вблизи несущей")
    ax[1].legend(fontsize=8, ncol=3, loc="upper center")
    ax[1].set_ylim(0.9, 1.16)
    fig.savefig(out / "fig_ach.png", dpi=170)
    plt.close(fig)

    fig, ax = new_fig(2, 3.1)
    ax[0].plot(f, np.unwrap(np.angle(ac["xd.v"])), "k--", lw=1.2, label="элемент задержки 5Δt")
    y = ac["y.v"]
    keep = np.abs(y) > 1e-6
    ax[0].plot(f[keep], np.unwrap(np.angle(y[keep])), color=VARIANTS[0]["color"], lw=1.2, label=VARIANTS[0]["name"])
    w = ac["yw.v"]
    keep = np.abs(w) > 1e-6
    ax[0].plot(f[keep], np.unwrap(np.angle(w[keep])), color=VARIANTS[1]["color"], lw=1.2, ls=":", label=VARIANTS[1]["name"])
    style_axes(ax[0], "f, МГц", "Φ, рад", "ФЧХ (развёрнутая) фильтров порядка 10 и элемента задержки")
    ax[0].set_xlim(0, 860)
    ax[0].legend(fontsize=8)
    for v in VARIANTS:
        d = ac[v["ac"] + ".v"] / ac[v["acref"] + ".v"]
        keep = np.abs(d) > 1e-6
        ax[1].plot(f[keep], np.angle(d[keep], deg=True), color=v["color"], lw=1.2, label=v["name"])
    ax[1].axvline(F0 / 1e6, color="#888888", ls=":", lw=1)
    ax[1].set_ylim(-100, -80)
    style_axes(ax[1], "f, МГц", "Δφ, град", "Разность фаз выхода фильтра и своего элемента задержки")
    ax[1].legend(fontsize=8)
    fig.savefig(out / "fig_fch.png", dpi=170)
    plt.close(fig)


def plot_signals(t, tr, m, tn, out):
    fig, ax = new_fig(2, 3.0)
    ax[0].plot(t * 1e6, tr["x0.Vt"], color="#2ca02c", lw=0.4)
    style_axes(ax[0], "t, мкс", "U, В", "Входной АМ-сигнал s(t): f0 = 430 МГц, Fм = 3 МГц, m = 1")
    win = (1.0e-6, 1.01e-6)
    k = (t >= win[0]) & (t <= win[1])
    kn = (tn >= win[0]) & (tn <= win[1])
    v = m["v"]
    ax[1].plot(t[k] * 1e9, tr["x0.Vt"][k], color="#2ca02c", lw=1, label="s(t)")
    ax[1].plot(t[k] * 1e9, tr[v["ref"] + ".Vt"][k], "k", lw=1.4, label="s'(t) = s(t − t0)")
    ax[1].plot(t[k] * 1e9, tr[v["out"] + ".Vt"][k], color=v["color"], ls="--", lw=1.4, label="ŝ(t)")
    ax[1].plot(tn[kn] * 1e9, m["sd"][kn], "ko", ms=3.5, label="s'(n)")
    ax[1].plot(tn[kn] * 1e9, m["sh"][kn], "s", color=v["color"], ms=3.5, label="ŝ(n)")
    style_axes(ax[1], "t, нс", "U, В", f"Осциллограммы: {v['name']}, t0 = 7Δt = 4,070 нс")
    ax[1].legend(fontsize=8, ncol=5, loc="upper center")
    ax[1].set_ylim(top=ax[1].get_ylim()[1] * 1.35)
    fig.savefig(out / "fig_signals.png", dpi=170)
    plt.close(fig)


def plot_no_delay(t, tr, tn, m, out):
    """Обоснование линии задержки: без неё z(t) = s(t) + jŝ(t) не аналитический."""
    v = m["v"]
    tw = (T_AN[0], T_AN[0] + 1 / FM)
    kn = (tn >= tw[0]) & (tn < tw[1])
    s_raw = np.interp(tn, t, tr["x0.Vt"])
    fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.9), constrained_layout=True)
    for a, re_part, title in (
        (ax[0], m["sd"], "с линией задержки: s'(n) + jŝ(n)"),
        (ax[1], s_raw, "без линии задержки: s(n) + jŝ(n)"),
    ):
        a.plot(re_part[kn], m["sh"][kn], ".", color=v["color"], ms=2.5)
        a.set_aspect("equal", "box")
        lim = 2.2
        a.set_xlim(-lim, lim)
        a.set_ylim(-lim, lim)
        style_axes(a, "Re z", "Im z", title)
    fig.savefig(out / "fig_no_delay.png", dpi=170)
    plt.close(fig)


def plot_hodograph(t, tr, m, tn, out):
    v = m["v"]
    k = (t >= T_AN[0]) & (t < T_AN[0] + HOD_FM_PERIODS / FM)
    kn = (tn >= T_AN[0]) & (tn < T_AN[1])
    sd, sh = tr[v["ref"] + ".Vt"], tr[v["out"] + ".Vt"]
    fig, ax = plt.subplots(1, 2, figsize=(8.2, 4.1), constrained_layout=True)
    ax[0].plot(sd[k], sh[k], color="#aaaaaa", lw=0.3, label="z(t)")
    ax[0].plot(m["z"][kn].real, m["z"][kn].imag, ".", color=v["color"], ms=4, label="z(n)")
    a = np.linspace(0, 2 * np.pi, 361)
    ax[0].plot(2 * np.cos(a), 2 * np.sin(a), "k:", lw=1, label="|z| = U0(1 + m)")
    ax[0].set_aspect("equal", "box")
    style_axes(ax[0], "s'(n), В", "ŝ(n), В", "Годограф z(n)")
    ax[0].legend(fontsize=8, loc="upper right", bbox_to_anchor=(1, 1.0))
    ax[1].plot(m["I"][kn], m["Q"][kn], ".", color=v["color"], ms=4)
    ax[1].set_xlim(-2.2, 2.2)
    ax[1].set_ylim(-2.2, 2.2)
    ax[1].set_aspect("equal", "box")
    style_axes(ax[1], "I(n), В", "Q(n), В", "Годограф в координатах I-Q")
    fig.savefig(out / "fig_hodograph.png", dpi=170)
    plt.close(fig)


def stem_db(ax, f, s, ref, color, label=None):
    y = np.maximum(db(s / ref), -100)
    ax.vlines(f, -100, y, color=color, lw=1.0, label=label)


def plot_spectra(sn, m, ian, fan, out):
    ss = spectrum(sn[ian])
    sz = spectrum(m["z"][ian])
    ref = sz.max()
    fig, ax = new_fig(3, 2.5)
    stem_db(ax[0], fan / 1e6, ss, ref, "#2ca02c")
    style_axes(ax[0], "f, МГц", "дБ", "Спектр исходного сигнала s(n): двусторонний")
    stem_db(ax[1], fan / 1e6, sz, ref, m["v"]["color"])
    style_axes(ax[1], "f, МГц", "дБ", f"Спектр аналитического сигнала z(n): подавление −f0 {-m['imgZ']:.1f} дБ".replace(".", ","))
    for a in ax:
        a.set_ylim(-100, 5)
        a.set_xlim(-860, 860)
    stem_db(ax[2], fan / 1e6, sz, ref, m["v"]["color"])
    ax[2].set_xlim(420, 440)
    style_axes(ax[2], "f, МГц", "дБ", "z(n) вблизи f0: несущая и боковые f0 ± 3 МГц (−6 дБ)")
    fig.savefig(out / "fig_spectra.png", dpi=170)
    plt.close(fig)


def plot_iq(m, tn, ian, fan, out):
    fig, ax = new_fig(3, 2.6)
    t0 = m["v"]["order"] / 2 * DT
    b = -1j * np.exp(-2j * np.pi * F0 * t0) * m["Uenv"]
    ax[0].plot(tn * 1e6, m["I"], color="#d62728", lw=1, label="I(n)")
    ax[0].plot(tn * 1e6, m["Q"], color="#1f77b4", lw=1, label="Q(n)")
    ax[0].plot(tn * 1e6, b.real, "k--", lw=1, label="I теор.")
    ax[0].plot(tn * 1e6, b.imag, "k:", lw=1, label="Q теор.")
    ax[0].set_xlim(0, tn[-1] * 1e6)
    style_axes(ax[0], "t, мкс", "U, В", "НЧ квадратурные составляющие I(n), Q(n)")
    ax[0].legend(fontsize=8, ncol=4, loc="upper right")
    si, sq = spectrum(m["I"][ian]), spectrum(m["Q"][ian])
    ref = max(si.max(), sq.max())
    stem_db(ax[1], fan / 1e6, si, ref, "#d62728", "|I(f)|")
    stem_db(ax[1], fan / 1e6 + 1.5, sq, ref, "#1f77b4", "|Q(f)| (смещён на 1,5 МГц для наглядности)")
    ax[1].set_ylim(-100, 5)
    ax[1].set_xlim(-860, 860)
    style_axes(ax[1], "f, МГц", "дБ", "Амплитудные спектры I и Q относительно постоянной составляющей")
    ax[1].legend(fontsize=8)
    stem_db(ax[2], fan / 1e6, si, ref, "#d62728")
    ax[2].set_xlim(-15, 15)
    ax[2].set_ylim(-100, 5)
    style_axes(ax[2], "f, МГц", "дБ", "Спектр I вблизи нуля: постоянная составляющая и ±3 МГц (−6 дБ)")
    fig.savefig(out / "fig_iq.png", dpi=170)
    plt.close(fig)


def plot_demod(m, tn, errs, out):
    fig, ax = new_fig(2, 3.0)
    ax[0].plot(tn * 1e6, m["Uenv"], "k--", lw=1.4, label="модулирующее U0(1 + m·sin 2πFм(t − t0))")
    ax[0].plot(tn * 1e6, m["U"], color=m["v"]["color"], lw=0.9, label="U(n)")
    ax[0].axvline(T_SKIP * 1e6, color="#888888", ls=":", lw=1)
    ax[0].set_xlim(0, tn[-1] * 1e6)
    style_axes(ax[0], "t, мкс", "U, В", f"Восстановленная огибающая, {m['v']['name']}")
    ax[0].legend(fontsize=8, loc="upper right")
    win = (tn >= 1e-6) & (tn <= 1.05e-6)
    for v in VARIANTS:
        ax[1].plot(tn[win] * 1e6, errs[v["key"]][win], color=v["color"], lw=0.9, label=v["name"])
    style_axes(ax[1], "t, мкс", "ΔU, В", "Ошибка U(n) − модулирующее (фрагмент)")
    ax[1].legend(fontsize=8, ncol=2)
    fig.savefig(out / "fig_demod.png", dpi=170)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dat", default="qucs", help="папка с .dat")
    ap.add_argument("--out", default="report/figures", help="папка для рисунков и результатов")
    a = ap.parse_args()
    dat, out = Path(a.dat), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})

    # --- AC (п. 3, 4) ---
    ac_names = ["acfrequency", "xd.v", "xd14.v", "xd18.v", "y.v", "yw.v", "y14.v", "y18.v"]
    ac = read_dat(dat / "LR2k5_N1_hilbert_AC.dat", ac_names)
    f = ac["acfrequency"]
    if2 = int(np.argmin(np.abs(f - F0)))
    plot_ac(ac, out)

    # теория для всех порядков (п. 4): выбор порядка
    orders = []
    for order in (10, 14, 18, 22):
        for ham in (False, True):
            h = coeffs(order, ham)
            k_f0 = abs(freq_resp(h, F0)[0])
            orders.append(dict(order=order, hamming=ham, K=k_f0, eps=(k_f0 - 1) * 100,
                               image=float(db(abs((1 - k_f0) / (1 + k_f0)))), taps=int(np.count_nonzero(h))))

    # --- TR (п. 5-8) ---
    tr_names = ["time", "x0.Vt", "xd.Vt", "xd14.Vt", "xd18.Vt", "y.Vt", "yw.Vt", "y14.Vt", "y18.Vt"]
    tr = read_dat(dat / "LR2k5_N1_hilbert_TR.dat", tr_names)
    t = tr["time"]
    s = tr["x0.Vt"]
    tn = np.arange(int(np.floor(t[-1] / DT + 1e-9)) + 1) * DT     # t_n = (n - 1) Δt
    ian = np.nonzero((tn >= T_AN[0] - DT / 2) & (tn < T_AN[1] - DT / 2))[0]
    fan = np.arange(-len(ian) // 2, len(ian) // 2) * FS / len(ian)
    sn = np.interp(tn, t, s)
    spec, freq = hilbert_ref(t, s)
    steady = t >= T_SKIP
    ok = tn >= T_SKIP

    def at(sp, f1):
        return sp[int(np.argmin(np.abs(fan - f1)))]

    res, errs, main_v = [], {}, None
    for v in VARIANTS:
        t0 = v["order"] / 2 * DT
        t0m, kt = measure_delay(tr[v["out"] + ".Vt"], spec, freq, steady, 20 * DT)
        sd = np.interp(tn, t, tr[v["ref"] + ".Vt"])
        sh = np.interp(tn, t, tr[v["out"] + ".Vt"])
        z = sd + 1j * sh
        wt = 2 * np.pi * F0 * tn                                   # ω0 (n - 1) Δt
        i_ = sd * np.cos(wt) + sh * np.sin(wt)
        q_ = -sd * np.sin(wt) + sh * np.cos(wt)
        u = np.sqrt(i_**2 + q_**2)
        uenv = U0 * (1 + MOD * np.sin(2 * np.pi * FM * (tn - t0)))
        sz = spectrum(z[ian])
        siq = spectrum((i_ + 1j * q_)[ian])
        yac, xac = ac[v["ac"] + ".v"], ac[v["acref"] + ".v"]
        ph = np.unwrap(np.angle(yac[if2 - 5:if2 + 6]))
        tg = -(ph[-1] - ph[0]) / (2 * np.pi * (f[if2 + 5] - f[if2 - 5]))
        r = dict(key=v["key"], name=v["name"], order=v["order"], t0_theory_ns=t0 * 1e9,
                 t0_ac_ns=tg * 1e9, t0_tr_ns=t0m * 1e9, K_ac=float(abs(yac[if2])), K_tr=float(kt),
                 dK_percent=float((abs(yac[if2]) - abs(xac[if2])) * 100),
                 dphi_ac=float(np.angle(yac[if2] / xac[if2], deg=True)),
                 dphi_tr=float(phase_at(sh[ian], sd[ian], tn[ian], F0)),
                 image_z_db=float(db(at(sz, -F0) / at(sz, F0))),
                 image_iq_db=float(db(at(siq, -FS / 2) / at(siq, 0))),
                 u_err_rms_percent=float(rms(u[ok] - uenv[ok]) / rms(uenv[ok]) * 100),
                 u_err_max=float(np.max(np.abs(u[ok] - uenv[ok]))),
                 q_max=float(np.max(np.abs(q_[ian]))), i_max=float(np.max(np.abs(i_[ian]))))
        res.append(r)
        errs[v["key"]] = u - uenv
        if v["key"] == MAIN:
            main_v = dict(v=v, sd=sd, sh=sh, z=z, I=i_, Q=q_, U=u, Uenv=uenv, imgZ=r["image_z_db"])

    sz_m = spectrum(main_v["z"][ian])
    si = spectrum(main_v["I"][ian])
    lines = {
        "fan_ref": [
            {"f_MHz": fm, "z_dB": float(db(at(sz_m, fm * 1e6) / sz_m.max()))}
            for fm in (-433, -430, -427, 427, 430, 433)
        ],
        "I_lines": [
            {"f_MHz": fm, "A": float(at(si, fm * 1e6))}
            for fm in (0, 3, -3, 860, -860)
        ],
    }

    plot_signals(t, tr, main_v, tn, out)
    plot_hodograph(t, tr, main_v, tn, out)
    plot_no_delay(t, tr, tn, main_v, out)
    plot_spectra(sn, main_v, ian, fan, out)
    plot_iq(main_v, tn, ian, fan, out)
    plot_demod(main_v, tn, errs, out)

    np.savetxt(out / "hilbert_iq.csv",
               np.column_stack([tn, sn, main_v["sd"], main_v["sh"], main_v["I"], main_v["Q"], main_v["U"], main_v["Uenv"]]),
               delimiter=";", header="t;s;s_d;s_h;I;Q;U;Uenv", comments="", fmt="%.6e")

    summary = dict(M=M, f0=F0, Fm=FM, fs=FS, dt_ns=DT * 1e9, f_at_table=float(f[if2]),
                   h10=coeffs(10, False).round(4).tolist(), h10w=coeffs(10, True).round(4).tolist(),
                   h14w=coeffs(14, True).round(4).tolist(), h18w=coeffs(18, True).round(4).tolist(),
                   ac_vs_theory_max=[float(np.max(np.abs(ac[v["ac"] + ".v"] - freq_resp(coeffs(v["order"], v["key"] != "o10"), f))))
                                     for v in VARIANTS],
                   orders=orders, results=res, lines=lines, tr_step_ps=float(np.median(np.diff(t)) * 1e12))
    (out / "results.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding="utf-8")

    sys.stdout.reconfigure(encoding="utf-8")
    print(f"fs = {FS / 1e6:.0f} МГц, Δt = {DT * 1e9:.4f} нс")
    for r in res:
        print(f"{r['name']:24s} t0 теор {r['t0_theory_ns']:.3f} AC {r['t0_ac_ns']:.3f} TR {r['t0_tr_ns']:.3f} нс | "
              f"|K| {r['K_ac']:.4f}/{r['K_tr']:.4f} | Δφ {r['dphi_ac']:.2f}/{r['dphi_tr']:.2f}° | "
              f"z(-f0) {r['image_z_db']:.1f} дБ | I/Q(fs/2) {r['image_iq_db']:.1f} дБ | СКО U {r['u_err_rms_percent']:.2f} % | max|Q| {r['q_max']:.3f}")


if __name__ == "__main__":
    main()
