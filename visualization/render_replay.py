"""Replay videos of the simulation, drawn from the data written by export_replay_data.m.

    python visualization/render_replay.py mission            # results/replay/mission.csv -> mission.mp4
    python visualization/render_replay.py outage             # results/replay/outage.csv  -> outage.mp4
    python visualization/render_replay.py mission --preview  # a few still frames only

Needs numpy, scipy, matplotlib and ffmpeg. Nothing is simulated here: attitude, rates, momentum and
the estimator bound are the logged signals of models/smallsat_adcs.slx. Only the playback speed is
chosen (time lapse, slower when the satellite turns fast), and it is printed on screen.
"""
import argparse
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial.transform import Rotation as Rot
from scipy.spatial.transform import Slerp

W, H, DPI, FPS = 1920, 1080, 120, 30
ARC = 180.0 / np.pi * 3600.0
PERIOD = 5801.0  # s, 600 km circular orbit

BG = "#080d1a"
PANEL = "#0d1426"
GRID = "#222c47"
TEXT = "#e9edf7"
MUTED = "#8b96b5"
ORANGE = "#ff7a59"
AMBER = "#ffc857"
GREEN = "#3ddc97"
BLUE = "#5aa9ff"
PURPLE = "#b28dff"
PINK = "#ff8ad4"
RED = "#ff5470"

for name in ("Inter", "Inter Display", "DejaVu Sans"):
    try:
        font_manager.findfont(name, fallback_to_default=False)
        plt.rcParams["font.family"] = name
        break
    except Exception:
        continue


# ----------------------------------------------------------------------------- data helpers
def load(path):
    return np.loadtxt(path, delimiter=",", skiprows=1)


def rot(q):
    """scalar-first quaternion rows -> scipy rotations (active), q = inertial to body"""
    q = np.asarray(q, float)
    return Rot.from_quat(np.column_stack([q[:, 1:4], q[:, 0]]))


def hms(s):
    s = int(round(s))
    return "%02d:%02d:%02d" % (s // 3600, (s // 60) % 60, s % 60)


def ease(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


class Frames:
    """pipes RGBA frames to ffmpeg, or keeps chosen frames as PNG"""

    def __init__(self, fig, out, preview_dir=None):
        self.fig, self.out, self.preview_dir = fig, out, preview_dir
        self.proc = None
        if preview_dir is None:
            cmd = ["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgba",
                   "-s", "%dx%d" % (W, H), "-r", str(FPS), "-i", "-", "-c:v", "libx264",
                   "-preset", "slow", "-crf", "17", "-pix_fmt", "yuv420p",
                   "-movflags", "+faststart", out]
            self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)

    def write(self, index, keep=False):
        self.fig.canvas.draw()
        if self.proc is not None:
            self.proc.stdin.write(np.asarray(self.fig.canvas.buffer_rgba()).tobytes())
        elif keep:
            self.fig.savefig(os.path.join(self.preview_dir, "frame_%04d.png" % index), dpi=DPI,
                             facecolor=BG)

    def close(self):
        if self.proc is not None:
            self.proc.stdin.close()
            self.proc.wait()


def style_axes(ax, ylabel, xlabel=None):
    ax.set_facecolor(PANEL)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(True, color=GRID, lw=0.6, alpha=0.8)
    ax.set_ylabel(ylabel, color=TEXT, fontsize=10)
    if xlabel:
        ax.set_xlabel(xlabel, color=MUTED, fontsize=9)


def background(fig, seed=3):
    ax = fig.add_axes([0, 0, 1, 1], zorder=-5)
    ax.set_facecolor(BG)
    r = np.random.RandomState(seed)
    n = 260
    ax.scatter(r.rand(n), r.rand(n), s=r.rand(n) * 2.2 + 0.2, c="white", alpha=0.35, lw=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    return ax


def overlay(fig):
    """full-figure layer for the intro and outro cards"""
    ov = fig.add_axes([0, 0, 1, 1], zorder=100)
    ov.axis("off")
    ov.set_xlim(0, 1)
    ov.set_ylim(0, 1)
    veil = ov.add_patch(plt.Rectangle((0, 0), 1, 1, color=BG, alpha=1.0, lw=0))
    return ov, veil


# ----------------------------------------------------------------------------- satellite model
def box_faces(lo, hi):
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    v = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                  [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
    idx = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [1, 2, 6, 5], [0, 3, 7, 4]]
    nrm = [[0, 0, -1], [0, 0, 1], [0, -1, 0], [0, 1, 0], [1, 0, 0], [-1, 0, 0]]
    return [v[i] for i in idx], np.array(nrm, float)


def build_satellite():
    """polygons (body frame), their normals and base colours"""
    polys, normals, colors = [], [], []
    f, n = box_faces((-0.5, -0.5, -0.7), (0.5, 0.5, 0.7))
    for p, nn in zip(f, n):
        polys.append(p)
        normals.append(nn)
        colors.append((0.80, 0.82, 0.88))
    # solar wings along +-x, plane normal y (faces the Sun in a dawn-dusk orbit)
    ncx, ncz = 5, 2
    for side in (1, -1):
        for i in range(ncx):
            for j in range(ncz):
                xa, xb = 0.62 + i * 0.44, 0.62 + (i + 1) * 0.44 - 0.03
                za, zb = -0.62 + j * 0.62, -0.62 + (j + 1) * 0.62 - 0.03
                p = np.array([[side * xa, 0, za], [side * xb, 0, za], [side * xb, 0, zb], [side * xa, 0, zb]])
                polys.append(p)
                normals.append(np.array([0.0, 1.0, 0.0]))
                colors.append((0.10, 0.22, 0.52))
        arm = np.array([[side * 0.5, -0.02, -0.04], [side * 0.62, -0.02, -0.04],
                        [side * 0.62, 0.02, 0.04], [side * 0.5, 0.02, 0.04]])
        polys.append(arm)
        normals.append(np.array([0.0, 1.0, 0.0]))
        colors.append((0.6, 0.6, 0.65))
    # payload aperture on +z (the nadir-facing side)
    a = np.linspace(0, 2 * np.pi, 28, endpoint=False)
    lens = np.column_stack([0.30 * np.cos(a), 0.30 * np.sin(a), np.full_like(a, 0.705)])
    polys.append(lens)
    normals.append(np.array([0.0, 0.0, 1.0]))
    colors.append((0.15, 0.85, 0.75))
    return polys, np.array(normals), np.array(colors)


def earth_cap(ax, zdepth):
    R = 9.0
    cz = zdepth - R
    th = np.linspace(0, np.radians(22), 40)
    ph = np.linspace(0, 2 * np.pi, 90)
    T, P = np.meshgrid(th, ph)
    X, Y, Z = R * np.sin(T) * np.cos(P), R * np.sin(T) * np.sin(P), cz + R * np.cos(T)
    shade = 0.5 + 0.5 * np.cos(T * 2.2)
    cols = np.zeros(T.shape + (4,))
    cols[..., 0] = 0.03 + 0.05 * shade
    cols[..., 1] = 0.12 + 0.30 * shade
    cols[..., 2] = 0.30 + 0.45 * shade
    cols[..., 3] = 1.0
    ax.plot_surface(X, Y, Z, facecolors=cols, shade=False, rstride=1, cstride=1, lw=0, antialiased=False, zorder=0)
    for lat in np.radians([4, 8, 12, 16, 20]):
        c = np.linspace(0, 2 * np.pi, 120)
        ax.plot(R * np.sin(lat) * np.cos(c), R * np.sin(lat) * np.sin(c), cz + R * np.cos(lat) + 0 * c,
                color=(0.5, 0.75, 1.0, 0.18), lw=0.6)
    for lon in np.radians(np.arange(0, 360, 30)):
        t = np.linspace(0, np.radians(22), 20)
        ax.plot(R * np.sin(t) * np.cos(lon), R * np.sin(t) * np.sin(lon), cz + R * np.cos(t),
                color=(0.5, 0.75, 1.0, 0.18), lw=0.6)


# ----------------------------------------------------------------------------- mission video
def render_mission(data_file, out, preview):
    d = load(data_file)
    t, q, qr, w, h, m, tw = d[:, 0], d[:, 1:5], d[:, 5:9], d[:, 9:12], d[:, 12:15], d[:, 15:18], d[:, 18:21]
    n = len(t)
    rq, rr = rot(q), rot(qr)
    slq, slr = Slerp(t, rq), Slerp(t, rr)
    rel = rr.inv() * rq
    err = np.degrees(rel.magnitude())
    rate = np.degrees(np.linalg.norm(w, axis=1))
    hn = np.linalg.norm(h, axis=1)
    dip = np.max(np.abs(m), axis=1)

    iE = np.flatnonzero(np.any(tw != 0, axis=1))
    tEng = t[iE[0]] if len(iE) else np.nan
    bad = np.flatnonzero((err >= 0.1) & (t >= tEng))
    tCap = t[bad[-1] + 1] if len(bad) and bad[-1] + 1 < n else tEng
    tail = t >= t[-1] - 2000
    mAfter = t >= tEng + 10
    stats = dict(tEng=tEng, tCap=tCap, mean=err[tail].mean(), mx=err[tail].max(),
                 dip=dip[mAfter].max(), h=hn[tail].max(), t0rate=rate[0])
    print("handover %.0f s, capture %.0f s after handover, pointing mean %.4f max %.4f deg, "
          "dipole peak %.1f A m^2" % (tEng, tCap - tEng, stats["mean"], stats["mx"], stats["dip"]))

    # time-lapse schedule: a frame every time the satellite has turned maxDeg, or at most dtMax apart
    fine = np.median(np.diff(t)) <= 0.5 + 1e-9      # data logged every 0.5 s: no interpolation needed in the close-up
    gstep = 0.5 if fine else 1.0
    tg = np.arange(t[0], t[-1] + 1e-9, gstep)
    rg = (slr(tg).inv() * slq(tg))
    step = (rg[1:] * rg[:-1].inv()).magnitude()
    step = np.degrees(step)
    dtMin, dtMax, target = 1.0, 40.0, 66.0   # seconds of video for the run itself
    # slow motion around the wheel capture: 0.5 s of flight per frame (x15) with 0.5 s data, else 1.5 s (x45)
    slowDt = 0.5 if fine else 1.5
    slowA, slowB = tEng - 10.0, tCap + (30.0 if fine else 120.0)

    def is_slow(ts):
        return slowA <= ts <= slowB

    def schedule(maxDeg):
        out_, acc, last = [tg[0]], 0.0, tg[0]
        for i in range(len(step)):
            acc += step[i]
            el = tg[i + 1] - last
            if is_slow(tg[i + 1]):
                emit = el >= slowDt
            else:
                emit = (acc >= maxDeg and el >= dtMin) or el >= dtMax
            if emit:
                out_.append(tg[i + 1])
                last, acc = tg[i + 1], 0.0
        if out_[-1] < t[-1]:
            out_.append(t[-1])
        return np.array(out_)

    lo, hi = 1.0, 60.0
    for _ in range(18):
        mid = np.sqrt(lo * hi)
        if len(schedule(mid)) > target * FPS:
            lo = mid
        else:
            hi = mid
    sched = schedule(hi)
    print("time lapse: a frame per %.1f deg of rotation (at most %.0f s apart)" % (hi, dtMax))
    intro, outro = int(2.4 * FPS), int(5.5 * FPS)
    nf = len(sched)
    iA, iB = np.searchsorted(sched, slowA), np.searchsorted(sched, slowB)
    print("slow motion from %.1f s to %.1f s of video (%d frames)" % ((intro + iA) / FPS, (intro + iB) / FPS, iB - iA))
    print("%d content frames (%.1f s), +%.1f s cards" % (nf, nf / FPS, (intro + outro) / FPS))
    warp = np.gradient(sched) * FPS
    warp = np.convolve(warp, np.ones(15) / 15, mode="same")

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI, facecolor=BG)
    background(fig)
    ax3 = fig.add_axes([-0.02, 0.20, 0.64, 0.68], projection="3d", computed_zorder=False)
    ax3.set_facecolor("none")
    ax3.set_axis_off()
    lim = 3.3
    ax3.set_xlim(-lim, lim)
    ax3.set_ylim(-lim, lim)
    ax3.set_zlim(-4.2, 2.4)
    ax3.set_box_aspect((2 * lim, 2 * lim, 6.6), zoom=1.2)
    ax3.view_init(elev=17, azim=-58)
    zE = -3.7
    earth_cap(ax3, zE)
    ax3.plot([0, 0], [0, 0], [0, zE], color=(1, 1, 1, 0.55), lw=1.2, ls=(0, (3, 3)), zorder=2)
    c = np.linspace(0, 2 * np.pi, 60)
    ax3.plot(0.35 * np.cos(c), 0.35 * np.sin(c), zE + 0 * c, color=GREEN, lw=1.6, alpha=0.9, zorder=2)
    ax3.plot([0], [0], [zE], marker="+", color=GREEN, ms=9, zorder=2)

    polys0, normals0, colors0 = build_satellite()
    coll = Poly3DCollection([p for p in polys0], linewidths=0.4, edgecolors=(0, 0, 0, 0.55), zorder=5)
    ax3.add_collection3d(coll)
    light = np.array([0.45, -0.35, 0.82])
    light /= np.linalg.norm(light)
    beam, = ax3.plot([0, 0], [0, 0], [0, zE], color=GREEN, lw=2.6, solid_capstyle="round", zorder=4)
    tip = ax3.scatter([0], [0], [zE], s=40, c=GREEN, zorder=4, depthshade=False)
    axl = [ax3.plot([0, 0], [0, 0], [0, 0], color=c_, lw=1.8, zorder=6)[0] for c_ in (RED, GREEN, BLUE)]
    axt = [ax3.text(0, 0, 0, s, color=c_, fontsize=9, fontweight="bold", zorder=7) for s, c_ in
           (("x", RED), ("y", GREEN), ("z", BLUE))]
    D = np.diag([1.0, -1.0, -1.0])  # show the orbit frame with zenith up: x along track, z up

    fig.text(0.035, 0.945, "Small-satellite attitude control", color=TEXT, fontsize=24, fontweight="bold",
             va="center")
    fig.text(0.035, 0.905, "tumble  →  B-dot detumble  →  reaction wheel nadir pointing  →  momentum dumping",
             color=MUTED, fontsize=11.5, va="center")
    fig.text(0.035, 0.872, "Orbit frame, Earth below. Dashed: nadir. Beam: where the payload axis (+z) actually points.",
             color=MUTED, fontsize=9.5)
    chip = fig.text(0.035, 0.035, "", fontsize=14, fontweight="bold", va="center",
                    bbox=dict(boxstyle="round,pad=0.45", fc=ORANGE, ec="none"), color="#101423")
    clock = fig.text(0.40, 0.030, "", color=TEXT, fontsize=13, va="center", family="DejaVu Sans Mono")
    lapse = fig.text(0.40, 0.072, "", color=MUTED, fontsize=9, va="center", linespacing=1.3)
    hud = []
    for i_, (lab, col_) in enumerate([("body rate", BLUE), ("pointing error", GREEN), ("wheel momentum", PURPLE)]):
        fig.text(0.04 + 0.17 * i_, 0.150, lab, color=MUTED, fontsize=10)
        hud.append(fig.text(0.04 + 0.17 * i_, 0.112, "", color=col_, fontsize=19, fontweight="bold", va="center"))

    # live plots
    specs = [("|ω| [deg/s]", BLUE, rate, (0.03, 10.0), True),
             ("pointing error [deg]", GREEN, np.maximum(err, 1e-4), (3e-3, 200.0), True),
             ("wheel momentum [N·m·s]", PURPLE, hn, (0, max(hn) * 1.18), False),
             ("torquer dipole [A·m²]", PINK, dip, (0, 33), False)]
    axs, lines, dots = [], [], []
    x = t / PERIOD
    # plotted lines: every 2 s of data, everything in and around the slow-motion window
    stride = max(1, int(round(2.0 / np.median(np.diff(t)))))
    selmask = np.zeros(len(t), bool)
    selmask[::stride] = True
    selmask |= (t >= slowA - 30.0) & (t <= slowB + 60.0)
    selmask[-1] = True
    sel = np.flatnonzero(selmask)
    tsel, xsel = t[sel], x[sel]
    ysel = [sp[2][sel] for sp in specs]
    top, hgt, gap = 0.905, 0.17, 0.055
    for i, (lab, col, y, yl, lg) in enumerate(specs):
        ax = fig.add_axes([0.635, top - (i + 1) * hgt - i * gap + gap * 0.0, 0.335, hgt])
        style_axes(ax, lab, "orbits" if i == 3 else None)
        if lg:
            ax.set_yscale("log")
        ax.set_xlim(0, x[-1])
        ax.set_ylim(*yl)
        ln, = ax.plot([], [], color=col, lw=1.6)
        dt, = ax.plot([], [], "o", color=col, ms=6)
        axs.append(ax)
        lines.append(ln)
        dots.append(dt)
    axs[1].axhline(0.1, color=AMBER, lw=1, ls="--")
    axs[1].text(x[-1] * 0.995, 0.115, "0.1°", color=AMBER, fontsize=8.5, ha="right", va="bottom")
    axs[3].axhline(30, color=RED, lw=1, ls="--")
    axs[3].text(x[-1] * 0.995, 28.5, "torquer limit 30", color=RED, fontsize=8.5, ha="right", va="top")
    evt = [[], [], [], []]
    fig.text(0.635, 0.935, "Logged signals of the Simulink model, estimate in the loop", color=MUTED, fontsize=10)

    ov, veil = overlay(fig)
    t1 = ov.text(0.5, 0.57, "Small-satellite ADCS", color=TEXT, fontsize=44, fontweight="bold", ha="center", va="center")
    t2 = ov.text(0.5, 0.47, "Simulink model, 12-state MEKF, B-dot detumble, wheel pointing, momentum dumping",
                 color=MUTED, fontsize=15, ha="center", va="center")
    t3 = ov.text(0.5, 0.40, "Generic Earth-observation smallsat, 600 km dawn-dusk Sun-synchronous orbit",
                 color=MUTED, fontsize=12, ha="center", va="center")
    outro_txt = [
        ov.text(0.5, 0.74, "Result of this run", color=MUTED, fontsize=15, ha="center"),
        ov.text(0.5, 0.64, "%.0f s from tip-off of %.1f°/s to wheel handover" % (tEng, stats["t0rate"]),
                color=TEXT, fontsize=22, ha="center"),
        ov.text(0.5, 0.55, "pointing locked to 0.1° within %.0f s of handover" % (tCap - tEng),
                color=TEXT, fontsize=22, ha="center"),
        ov.text(0.5, 0.46, "pointing error %.4f° mean, %.4f° max (last 2000 s)" % (stats["mean"], stats["mx"]),
                color=GREEN, fontsize=22, ha="center", fontweight="bold"),
        ov.text(0.5, 0.37, "torquer dipole peak %.0f A·m² of 30 after handover" % stats["dip"],
                color=TEXT, fontsize=22, ha="center"),
        ov.text(0.5, 0.22, "Open source, MATLAB / Simulink  |  all numbers are from the simulation, no hardware",
                color=MUTED, fontsize=13, ha="center"),
    ]
    for tx in [t1, t2, t3] + outro_txt:
        tx.set_alpha(0)

    def draw_state(ts, k, wp):
        R = (slr(ts).inv() * slq(ts))
        # satellite pose: body to orbit frame = R(ref)^T R(q)  ->  slr.inv * slq already is (body in ref)
        Rm = D @ rel_matrix(ts)
        pv = []
        for p in polys0:
            pv.append((Rm @ p.T).T)
        nrm = (Rm @ normals0.T).T
        sh = np.clip(nrm @ light, -1, 1) * 0.5 + 0.5
        cols = []
        for c_, s_, kind in zip(colors0, sh, range(len(colors0))):
            f = 0.35 + 0.75 * s_
            cols.append(np.clip(np.array(c_) * f, 0, 1).tolist() + [1.0])
        coll.set_verts(pv)
        coll.set_facecolor(cols)
        zb = Rm @ np.array([0, 0, 1.0])
        L = 3.7
        e_now = np.interp(ts, t, err)
        gcol = np.array(matplotlib.colors.to_rgb(GREEN))
        rcol = np.array(matplotlib.colors.to_rgb(ORANGE))
        a = float(np.clip(np.log10(max(e_now, 0.1) / 0.1) / 2.0, 0, 1))
        col = tuple(gcol * (1 - a) + rcol * a)
        beam.set_data_3d([0, zb[0] * L], [0, zb[1] * L], [0, zb[2] * L])
        beam.set_color(col)
        tip._offsets3d = ([zb[0] * L], [zb[1] * L], [zb[2] * L])
        tip.set_color([col])
        for i_, ev in enumerate(np.eye(3)):
            v = Rm @ ev * 1.45
            axl[i_].set_data_3d([0, v[0]], [0, v[1]], [0, v[2]])
            axt[i_].set_position_3d((v[0] * 1.08, v[1] * 1.08, v[2] * 1.08))
        # plots
        kk = np.searchsorted(tsel, ts, side="right")
        for ln, dt, ys_, (lab, colr, y, yl, lg) in zip(lines, dots, ysel, specs):
            ln.set_data(xsel[:kk], ys_[:kk])
            dt.set_data([ts / PERIOD], [np.interp(ts, t, y)])
        for i_ in range(4):
            if ts >= tEng and len(evt[i_]) == 0:
                evt[i_].append(axs[i_].axvline(tEng / PERIOD, color=AMBER, lw=1, ls=":"))
                if i_ == 0:
                    axs[i_].text(tEng / PERIOD, yl_top(axs[i_]), " wheels take over", color=AMBER, fontsize=8.5, va="top")
            if ts >= tCap and len(evt[i_]) == 1:
                evt[i_].append(axs[i_].axvline(tCap / PERIOD, color=GREEN, lw=1, ls=":"))
                if i_ == 0:
                    axs[i_].text(tCap / PERIOD, yl_top(axs[i_]) * 0.5, " pointing locked", color=GREEN, fontsize=8.5, va="top")
        if ts < tEng:
            chip.set_text("  DETUMBLE  magnetorquers, B-dot  ")
            chip.get_bbox_patch().set_facecolor(ORANGE)
        elif ts < tCap:
            chip.set_text("  CAPTURE  reaction wheels  ")
            chip.get_bbox_patch().set_facecolor(AMBER)
        else:
            chip.set_text("  NADIR POINTING  wheels + momentum dumping  ")
            chip.get_bbox_patch().set_facecolor(GREEN)
        clock.set_text("T+ " + hms(ts))
        if is_slow(ts):
            lapse.set_text("slow motion x%d\n(%s)" % (int(round(wp)), "logged every 0.5 s" if fine else "interpolated between 10 s log samples"))
            lapse.set_color(AMBER)
        else:
            lapse.set_text("time lapse x%d" % int(round(wp)))
            lapse.set_color(MUTED)
        r_ = np.interp(ts, t, rate)
        e_ = np.interp(ts, t, err)
        hud[0].set_text("%.2f \u00b0/s" % r_)
        hud[1].set_text(("%.2f\u00b0" % e_) if e_ >= 1 else ("%.4f\u00b0" % e_))
        hud[2].set_text("%.3f N\u00b7m\u00b7s" % np.interp(ts, t, hn))

    def rel_matrix(ts):
        r = slr(ts).inv() * slq(ts)   # body vectors -> reference-frame vectors
        return r.as_matrix()

    def yl_top(ax):
        return ax.get_ylim()[1]

    total = intro + nf + outro
    marks = {0, intro - 1, intro + nf // 4, intro + nf // 2, intro + int(nf * 0.62), intro + int(nf * 0.8),
             total - outro - 1, total - 20}
    pdir = None
    if preview:
        pdir = os.path.join(os.path.dirname(out), "preview")
        os.makedirs(pdir, exist_ok=True)
    fr = Frames(fig, out, pdir)
    draw_state(sched[0], 0, warp[0])
    for f in range(total):
        if f < intro:
            k = 0
            veil.set_alpha(1.0 - float(ease((f - 1.5 * FPS) / (0.8 * FPS))))
            for tx in (t1, t2, t3):
                tx.set_alpha(float(ease(f / (0.5 * FPS))) * float(1 - ease((f - 1.5 * FPS) / (0.6 * FPS))))
            draw_state(sched[0], 0, warp[0])
        elif f < intro + nf:
            k = f - intro
            veil.set_alpha(0)
            for tx in (t1, t2, t3):
                tx.set_alpha(0)
            draw_state(sched[k], k, warp[k])
        else:
            g = (f - intro - nf) / (0.9 * FPS)
            veil.set_alpha(0.95 * float(ease(g)))
            for i_, tx in enumerate(outro_txt):
                tx.set_alpha(float(ease(g - 0.15 * i_)))
            draw_state(sched[-1], nf - 1, warp[-1])
        if preview and f not in marks:
            continue
        fr.write(f, keep=f in marks)
        if f % 150 == 0:
            print("frame %d / %d" % (f, total), flush=True)
    fr.close()


# ----------------------------------------------------------------------------- outage video
def render_outage(data_file, out, preview, t_out0=14000.0, t_out1=15800.0):
    d = load(data_file)
    t = d[:, 0]
    q, qr, qh, sg = d[:, 1:5], d[:, 5:9], d[:, 9:13], d[:, 13:16]
    rq, rr, rh = rot(q), rot(qr), rot(qh)
    # knowledge error as in run_outage: rotation from the estimate to the truth, total angle
    e = np.degrees((rh.inv() * rq).magnitude()) * 3600.0
    bound = 3.0 * np.linalg.norm(sg, axis=1) * ARC
    pt = np.degrees((rr.inv() * rq).magnitude())
    ta, tb = t_out0 - 1200.0, t_out1 + 1200.0
    tg = np.arange(ta, tb + 1e-9, 5.0)
    f_ = lambda y: np.interp(tg, t, y)
    eG, bG, pG = f_(e), f_(bound), f_(pt)
    during = (tg >= t_out0) & (tg < t_out1)
    inside = float(np.mean(e[(t >= t_out0) & (t < t_out1)] < bound[(t >= t_out0) & (t < t_out1)]))
    peak = float(e[(t >= t_out0) & (t < t_out1)].max())
    peakB = float(bound[(t >= t_out0) & (t < t_out1)].max())
    late = np.flatnonzero((t >= t_out1) & (pt >= 0.03))
    recover = (t[late[-1]] + 10.0 - t_out1) if len(late) else 0.0
    ptDuring = float(pt[(t >= t_out0) & (t < t_out1)].max())
    print("outage: error peak %.0f arcsec, 3 sigma peak %.0f arcsec, inside %.1f %%, pointing peak %.3f deg, "
          "back under 0.03 deg %.0f s after the fixes return" % (peak, peakB, 100 * inside, ptDuring, recover))

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI, facecolor=BG)
    background(fig, seed=8)
    fig.text(0.04, 0.945, "Star tracker lost for 30 minutes", color=TEXT, fontsize=26, fontweight="bold", va="center")
    fig.text(0.04, 0.905, "The estimator runs on the gyro alone. Its own uncertainty tells how wrong the attitude can be.",
             color=MUTED, fontsize=12, va="center")
    a1 = fig.add_axes([0.065, 0.36, 0.60, 0.50])
    a2 = fig.add_axes([0.065, 0.095, 0.60, 0.20])
    style_axes(a1, "attitude knowledge error [arcsec]")
    style_axes(a2, "pointing error [deg]", "time [min from outage start]")
    xm = (tg - t_out0) / 60.0
    fig.text(0.065, 0.012, "Logged every 10 s, drawn with linear interpolation.", color=MUTED, fontsize=9)
    for a in (a1, a2):
        a.set_xlim(xm[0], xm[-1])
        a.axvspan(0, (t_out1 - t_out0) / 60.0, color=RED, alpha=0.09, lw=0)
    a1.set_yscale("log")
    a1.set_ylim(2, 6000)
    a2.set_yscale("log")
    a2.set_ylim(1e-3, 1.0)
    a2.axhline(0.03, color=AMBER, lw=1, ls="--")
    a2.text(xm[-1], 0.034, "0.03° ", color=AMBER, fontsize=8.5, ha="right", va="bottom")
    a1.text(15, 2.6, "STAR TRACKER BLIND", color=RED, fontsize=12, ha="center", fontweight="bold")
    l_b, = a1.plot([], [], color=PURPLE, lw=2.2, label="filter 3σ bound")
    fill = [None]
    l_e, = a1.plot([], [], color=GREEN, lw=1.8, label="true error")
    d_e, = a1.plot([], [], "o", color=GREEN, ms=7)
    d_b, = a1.plot([], [], "o", color=PURPLE, ms=7)
    leg = a1.legend(loc="upper left", frameon=False, labelcolor=TEXT, fontsize=11)
    l_p, = a2.plot([], [], color=BLUE, lw=1.6)
    d_p, = a2.plot([], [], "o", color=BLUE, ms=6)

    cx = 0.72
    chip = fig.text(cx, 0.80, "", fontsize=15, fontweight="bold", va="center", color="#101423",
                    bbox=dict(boxstyle="round,pad=0.5", fc=GREEN, ec="none"))
    cards = []
    for i, (lab, col) in enumerate([("true knowledge error", GREEN), ("filter 3σ bound", PURPLE),
                                    ("error inside the bound", TEXT)]):
        y0 = 0.64 - i * 0.14
        fig.text(cx, y0 + 0.045, lab, color=MUTED, fontsize=11)
        cards.append(fig.text(cx, y0, "", color=col, fontsize=30, fontweight="bold", va="center"))
    clock = fig.text(cx, 0.865, "", color=TEXT, fontsize=14, family="DejaVu Sans Mono", va="center")
    verdict = ("The filter does not hide this: the true error stays\ninside its bound." if inside >= 0.95
               else "In this run the true error leaves the bound in\n%.0f %% of samples." % (100 * (1 - inside)))
    fig.text(cx, 0.19, "Why it grows: gyro angle random walk, 0.003°/√s.\nError grows like √t, so the pointing law cannot be\nheld to 0.03° without star fixes after about 30 s.\n" + verdict,
             color=MUTED, fontsize=11, va="top", linespacing=1.5)

    ov, veil = overlay(fig)
    outro = [ov.text(0.5, 0.66, "Star tracker outage of 1800 s", color=MUTED, fontsize=16, ha="center"),
             ov.text(0.5, 0.56, "true error peaked at %.0f arcsec, bound %.0f arcsec" % (peak, peakB), color=TEXT, fontsize=24, ha="center"),
             ov.text(0.5, 0.47, "error inside the 3σ bound in %.1f %% of samples" % (100 * inside), color=GREEN, fontsize=24, fontweight="bold", ha="center"),
             ov.text(0.5, 0.38, ("pointing back under 0.03° within %.0f s of the fixes returning" % max(recover, 10.0)) if recover <= 60 else ("pointing under 0.03° again %.0f s after the fixes return" % recover), color=TEXT, fontsize=24, ha="center")]
    for tx in outro:
        tx.set_alpha(0)

    n = len(tg)
    stride = 1
    idx = list(range(0, n, stride))
    intro, outro_n = int(1.0 * FPS), int(5.0 * FPS)
    total = intro + len(idx) + outro_n
    pdir = None
    if preview:
        pdir = os.path.join(os.path.dirname(out), "preview")
        os.makedirs(pdir, exist_ok=True)
    marks = {intro + len(idx) // 4, intro + int(len(idx) * 0.45), intro + int(len(idx) * 0.62), intro + len(idx) - 2,
             total - 10}
    fr = Frames(fig, out, pdir)

    def draw(k):
        kk = k + 1
        l_b.set_data(xm[:kk], bG[:kk])
        l_e.set_data(xm[:kk], np.maximum(eG[:kk], 1e-3))
        d_e.set_data([xm[k]], [eG[k]])
        d_b.set_data([xm[k]], [bG[k]])
        l_p.set_data(xm[:kk], np.maximum(pG[:kk], 1e-4))
        d_p.set_data([xm[k]], [pG[k]])
        if fill[0] is not None:
            fill[0].remove()
        fill[0] = a1.fill_between(xm[:kk], 1e-3, bG[:kk], color=PURPLE, alpha=0.10, lw=0)
        blind = t_out0 <= tg[k] < t_out1
        chip.set_text("  STAR TRACKER  BLIND  " if blind else "  STAR TRACKER  LOCKED  ")
        chip.get_bbox_patch().set_facecolor(RED if blind else GREEN)
        clock.set_text("outage %s / 00:30:00" % hms(min(max(tg[k] - t_out0, 0), t_out1 - t_out0)))
        cards[0].set_text("%.1f arcsec" % eG[k])
        cards[1].set_text("%.1f arcsec" % bG[k])
        cards[2].set_text("yes" if eG[k] < bG[k] else "NO")

    for f in range(total):
        if f < intro:
            veil.set_alpha(1 - float(ease(f / (0.8 * FPS))))
            draw(0)
        elif f < intro + len(idx):
            veil.set_alpha(0)
            draw(idx[f - intro])
        else:
            g = (f - intro - len(idx)) / (0.9 * FPS)
            veil.set_alpha(0.95 * float(ease(g)))
            for i_, tx in enumerate(outro):
                tx.set_alpha(float(ease(g - 0.2 * i_)))
            draw(idx[-1])
        if preview and f not in marks:
            continue
        fr.write(f, keep=f in marks)
    fr.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=["mission", "outage"])
    ap.add_argument("--data", default=os.path.join(os.path.dirname(__file__), "..", "results", "replay"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--preview", action="store_true", help="write a few PNG frames instead of the video")
    a = ap.parse_args()
    src = os.path.join(a.data, a.which + ".csv")
    dst = a.out or os.path.join(a.data, a.which + ".mp4")
    (render_mission if a.which == "mission" else render_outage)(src, dst, a.preview)
    print("done", dst)
