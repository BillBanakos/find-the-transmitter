"""
"Find the transmitter" - RADAR edition, robust detection for CW *or* pulsed sources.
Source: phone Wi-Fi hotspot (PHONE_MODE = True) or spectrum analyzer tone (PHONE_MODE = False).
PlutoSDR+ (2r2t) + 2 monopoles. See README.md for setup and the legal notes on transmitting.

Detection runs two detectors on every buffer and keeps the better one:
  NB  (narrowband) : band-pass around the strongest spectral line   -> best for CW / long pulses
  TD  (time-domain): gate on raw envelope power above the noise     -> best for short pulses
Phase difference RX2 vs RX1 is taken only from the detected samples.

Right-hand panels: spectrum, RX envelope (you SEE the pulses), phasors.
Console prints a diagnostic line ~1x per second - send it if detection fails.

Anti-lag : 1 kernel buffer (no stale data), blitted drawing, heavy panels every 3rd frame.
Anti-jitter: frames with low RX1/RX2 coherence are rejected, then a 5-point median + light smoothing.
Tune MED_LEN / ALPHA / MIN_COH in the config if you want it steadier or snappier.

Keys: g = auto gain   c = calibrate at 0 deg   f = flip   m = Bartlett/MUSIC   r = clear track   q = quit
Order: source ON -> g -> volunteer at 0 deg -> c -> walk

Run:  python pluto_aoa_radar.py          (hardware)
      python pluto_aoa_radar.py --sim    (no hardware, simulated transmitter)
"""
from collections import deque
import time
import sys
import numpy as np

if "--sim" in sys.argv:          # no hardware: simulated Pluto + walking transmitter
    import sim_pluto as adi
else:
    import adi
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ---------------- config ----------------
URI = "ip:192.168.2.1"
PHONE_MODE = True               # True: phone Wi-Fi hotspot beacons (legal, certified TX). False: SA tone
WIFI_CH = 6                     # 2.4 GHz channel of the hotspot (find it with a Wi-Fi analyzer app / the SA)
if PHONE_MODE:
    F_SIG = 2407e6 + 5e6 * WIFI_CH  # channel centre
    FC = F_SIG                      # beacon is ~20 MHz wide, we listen to the middle 2 MHz
    N = 2**18                       # 131 ms > 102.4 ms beacon interval -> >=1 beacon per buffer
    GAIN0 = 10                      # phone is ~35 dB stronger than the SA
else:
    F_SIG = 2.478e9                 # SA tone: inside 2400-2483.5 MHz SRD band, above Wi-Fi ch 11
    FC = F_SIG - 200e3              # tone lands at +200 kHz, away from DC/LO leakage
    N = 2**15                       # 16.4 ms per buffer
    GAIN0 = 40
FS = 2.0e6
LAM = 3e8 / F_SIG
D = 0.0625                      # MEASURED element spacing [m] (lambda/2 at 2.478 GHz would be 6.05 cm)
FILT_BW = 100e3                 # NB detector half-bandwidth [Hz]
MIN_SNR_DB = 8                  # below this -> NO SIGNAL
MIN_SAMPLES = 8                 # min detected samples per buffer
MIN_COH = 0.85                  # phase coherence RX1/RX2 (1 = clean single source); below -> reject
MED_LEN = 5                     # median filter length (kills random jumps)
ALPHA = 0.5                     # extra smoothing after median (0 = none, 0.9 = heavy)
SLOW_EVERY = 3                  # update spectrum/envelope panels every Nth frame
FULL_SCALE = 2048               # 12-bit ADC
TRAIL = 60
DC_GUARD = int(20e3 / FS * N)
ENV_DEC = max(1, N // 2048)     # envelope plot decimation

GREEN, CYAN, MAG, AMBER, RED = "#39ff14", "#00e5ff", "#ff2fd6", "#ffb000", "#ff3030"

# ---------------- radio ----------------
sdr = adi.ad9361(uri=URI)
sdr.sample_rate = int(FS)
sdr.rx_lo = int(FC)
sdr.rx_rf_bandwidth = int(FS)
sdr.rx_enabled_channels = [0, 1]
sdr.rx_buffer_size = N
try:
    sdr._rxadc.set_kernel_buffers_count(1)       # no queued (stale) buffers -> no lag
except Exception as e:
    print("could not set kernel buffers:", e)
sdr.gain_control_mode_chan0 = "manual"
sdr.gain_control_mode_chan1 = "manual"
sdr.tx_hardwaregain_chan0 = -89
sdr.tx_hardwaregain_chan1 = -89


def set_gain(g):
    g = int(np.clip(g, 0, 71))
    sdr.rx_hardwaregain_chan0 = g
    sdr.rx_hardwaregain_chan1 = g
    return g


state = {"cal": 0.0, "sign": 1.0, "music": True, "run": True,
         "do_cal": False, "do_gain": False, "gain": set_gain(GAIN0)}

k_d = 2 * np.pi * D / LAM
th = np.radians(np.linspace(-90, 90, 361))
freqs = np.fft.fftfreq(N, 1 / FS)
win = np.hanning(N)
t_ms = np.arange(N // ENV_DEC) * ENV_DEC / FS * 1e3


def read_raw():
    x0, x1 = (np.asarray(x, dtype=np.complex128) for x in sdr.rx())
    return x0 - x0.mean(), x1 - x1.mean()          # remove DC / LO leakage


def auto_gain():
    """Put the strongest samples at ~25-60 % of full scale, then lock gain."""
    g, pk = state["gain"], 0.0
    for _ in range(15):
        read_raw()                                  # flush after gain change
        x0, x1 = read_raw()
        pk = max(np.abs(x0).max(), np.abs(x1).max()) / FULL_SCALE
        if pk > 0.6:
            g = set_gain(g - 6)
        elif pk < 0.25 and g < 71:
            g = set_gain(g + 3)
        else:
            break
    state["gain"] = g
    print(f"[gain] locked at {g} dB, peak {pk*100:.0f} % FS -> now press c at 0 deg")


def detect_nb(F0, F1):
    mag = np.abs(F0) + np.abs(F1)
    mag[:DC_GUARD] = 0
    mag[-DC_GUARD:] = 0
    k = int(np.argmax(mag))
    mask = np.abs(freqs - freqs[k]) < FILT_BW
    y0, y1 = np.fft.ifft(F0 * mask), np.fft.ifft(F1 * mask)
    p = np.abs(y0) ** 2 + np.abs(y1) ** 2
    out = ~mask
    out[:DC_GUARD] = False
    out[-DC_GUARD:] = False
    sigma2N = (np.median(np.abs(F0[out]) ** 2) + np.median(np.abs(F1[out]) ** 2)) / np.log(2)
    noise = sigma2N * mask.sum() / N**2
    on = p > max(0.25 * p.max(), 6 * noise)
    snr = 10 * np.log10(p[on].mean() / noise) if on.sum() >= MIN_SAMPLES else 0.0
    return snr, on, y0[on], y1[on], freqs[k] / 1e3


def detect_td(x0, x1, p):
    noise = np.percentile(p, 20) / -np.log(0.8)    # mean of exponential noise power
    on = p > max(0.25 * p.max(), 10 * noise)
    snr = 10 * np.log10(p[on].mean() / noise) if on.sum() >= MIN_SAMPLES else 0.0
    return snr, on, x0[on], x1[on], noise


def acquire(slow):
    x0, x1 = read_raw()
    F0, F1 = np.fft.fft(x0), np.fft.fft(x1)
    p = np.abs(x0) ** 2 + np.abs(x1) ** 2
    pk = np.sqrt(p.max()) / FULL_SCALE
    snr_nb, on_nb, a0, a1, f_pk = detect_nb(F0, F1)
    snr_td, on_td, b0, b1, noise_td = detect_td(x0, x1, p)
    if snr_nb >= snr_td:
        mode, snr, on, s0, s1 = "NB", snr_nb, on_nb, a0, a1
    else:
        mode, snr, on, s0, s1 = "TD", snr_td, on_td, b0, b1
    corr = np.sum(s1 * np.conj(s0))
    coh = np.abs(corr) / np.sqrt(np.sum(np.abs(s0) ** 2) * np.sum(np.abs(s1) ** 2) + 1e-20)
    d = dict(S=np.vstack([s0, s1]), corr=corr, coh=coh, snr=snr, snr_nb=snr_nb, snr_td=snr_td,
             mode=mode, duty=on.mean(), n_on=int(on.sum()), f_pk=f_pk, pk=pk)
    if slow:                                        # display-only work, not every frame
        W0, W1 = np.fft.fft(x0 * win), np.fft.fft(x1 * win)
        d["sp0"] = 20 * np.log10(np.abs(np.fft.fftshift(W0)) + 1e-9)
        d["sp1"] = 20 * np.log10(np.abs(np.fft.fftshift(W1)) + 1e-9)
        d["env_db"] = 10 * np.log10(p[: N // ENV_DEC * ENV_DEC].reshape(-1, ENV_DEC).max(1)
                                    / FULL_SCALE**2 + 1e-12)
        d["thr_db"] = 10 * np.log10(10 * noise_td / FULL_SCALE**2 + 1e-12)
    return d


def on_key(ev):
    flags = {"c": "do_cal", "g": "do_gain"}
    if ev.key in flags:
        state[flags[ev.key]] = True
    elif ev.key == "f":
        state["sign"] *= -1
        hist.clear()
    elif ev.key == "m":
        state["music"] = not state["music"]
    elif ev.key == "r":
        trail.clear()
    elif ev.key == "q":
        state["run"] = False


# ---------------- figure ----------------
plt.style.use("dark_background")
fig = plt.figure(figsize=(14, 8), facecolor="black")
fig.canvas.manager.set_window_title("IEEE SB Patras - SDR direction finding")
gs = GridSpec(3, 3, figure=fig, hspace=0.55, wspace=0.25, top=0.9, bottom=0.08)

axr = fig.add_subplot(gs[0:3, 0:2], projection="polar")
axr.set_facecolor("#001400")
axr.set_thetamin(-90); axr.set_thetamax(90)
axr.set_theta_zero_location("N"); axr.set_theta_direction(-1)
axr.set_rlim(0, 1.05); axr.set_yticklabels([])
axr.set_thetagrids(range(-90, 91, 15), color=GREEN, fontsize=10)
axr.grid(color=GREEN, alpha=0.25)
for r in (0.25, 0.5, 0.75, 1.0):
    axr.plot(th, np.full_like(th, r), color=GREEN, alpha=0.15, lw=0.8)
beam_fill = axr.fill(th, np.zeros_like(th), color=GREEN, alpha=0.18)[0]
beam_line, = axr.plot(th, np.zeros_like(th), color=GREEN, lw=2)
sweep, = axr.plot([0, 0], [0, 1.05], color=GREEN, lw=1.5, alpha=0.35)
glow, = axr.plot([0, 0], [0, 1], color=MAG, lw=12, alpha=0.25, solid_capstyle="round")
needle, = axr.plot([0, 0], [0, 1], color=MAG, lw=3.5, solid_capstyle="round")
trail = deque(maxlen=TRAIL)
trail_sc = axr.scatter([0], [0], s=0, c=[0], cmap="cool", vmin=0, vmax=1, edgecolors="none")
big = fig.text(0.36, 0.17, "", ha="center", fontsize=50, color=MAG, family="monospace", weight="bold")
info = fig.text(0.36, 0.10, "", ha="center", fontsize=11, color=GREEN, family="monospace")
fig.text(0.36, 0.94, "WHERE IS THE TRANSMITTER?", ha="center", fontsize=22, color=GREEN,
         family="monospace", weight="bold")

axs = fig.add_subplot(gs[0, 2])
axs.set_title("Spectrum (RX1 / RX2)", color=CYAN, fontsize=10)
SDEC = max(1, N // 2048)
fd = np.fft.fftshift(freqs)[: N // SDEC * SDEC].reshape(-1, SDEC).mean(1) / 1e3
sp0, = axs.plot(fd, np.full_like(fd, -100), color=CYAN, lw=1)
sp1, = axs.plot(fd, np.full_like(fd, -100), color=AMBER, lw=1, alpha=0.8)
axs.set_xlim(fd[0], fd[-1]); axs.set_ylim(-100, 5)
axs.set_xlabel("offset [kHz]", fontsize=8); axs.set_ylabel("dBc", fontsize=8)
axs.grid(alpha=0.2); axs.tick_params(labelsize=8)

axe = fig.add_subplot(gs[1, 2])
axe.set_title("RX envelope (pulses)", color=CYAN, fontsize=10)
env_line, = axe.plot(t_ms, np.full_like(t_ms, -100), color=GREEN, lw=1)
thr_line = axe.axhline(-60, color=RED, lw=1, ls="--", alpha=0.8)
axe.set_xlim(0, t_ms[-1]); axe.set_ylim(-90, 5)
axe.set_xlabel("time [ms]", fontsize=8); axe.set_ylabel("dBFS", fontsize=8)
axe.grid(alpha=0.2); axe.tick_params(labelsize=8)

axp = fig.add_subplot(gs[2, 2])
axp.set_title("Phasors: RX1 (ref) vs RX2", color=CYAN, fontsize=10)
c = np.linspace(0, 2 * np.pi, 200)
axp.plot(np.cos(c), np.sin(c), color="gray", alpha=0.4)
axp.plot([0, 1], [0, 0], color=CYAN, lw=3)
ph1, = axp.plot([0, 1], [0, 0], color=AMBER, lw=3)
ph_txt = axp.text(1.6, 0, "", ha="left", va="center", color="white", fontsize=11, family="monospace")
axp.set_xlim(-1.2, 3.6); axp.set_ylim(-1.2, 1.2); axp.set_aspect("equal"); axp.axis("off")

fig.canvas.mpl_connect("key_press_event", on_key)

# ---- blitting: static background is drawn once, only moving parts are redrawn each frame ----
dyn = [beam_fill, beam_line, sweep, glow, needle, trail_sc, big, info,
       sp0, sp1, env_line, thr_line, ph1, ph_txt]
for art in dyn:
    art.set_animated(True)
bg = {"img": None}


def on_draw(_):
    bg["img"] = fig.canvas.copy_from_bbox(fig.bbox)


fig.canvas.mpl_connect("draw_event", on_draw)
plt.show(block=False)
fig.canvas.draw()
fig.canvas.flush_events()


def render():
    if bg["img"] is None:
        fig.canvas.draw()
    fig.canvas.restore_region(bg["img"])
    for art in dyn:
        fig.draw_artist(art)
    fig.canvas.blit(fig.bbox)
    fig.canvas.flush_events()


# ---------------- loop ----------------
hist = deque(maxlen=MED_LEN)            # recent raw bearings for the median filter
est, sweep_ph, dphi, t_log, frame, fps = None, 0.0, 0.0, 0.0, 0, 0.0
t_prev = time.time()
while state["run"] and plt.fignum_exists(fig.number):
    if state["do_gain"]:
        state["do_gain"] = False
        auto_gain()
        hist.clear()
    slow = frame % SLOW_EVERY == 0
    a_ = acquire(slow)
    frame += 1

    sweep_ph = (sweep_ph + 0.12) % np.pi
    sweep.set_data([sweep_ph - np.pi / 2] * 2, [0, 1.05])
    if slow:
        top = max(a_["sp0"].max(), a_["sp1"].max())
        sp0.set_ydata((a_["sp0"] - top)[: N // SDEC * SDEC].reshape(-1, SDEC).max(1))
        sp1.set_ydata((a_["sp1"] - top)[: N // SDEC * SDEC].reshape(-1, SDEC).max(1))
        env_line.set_ydata(a_["env_db"])
        thr_line.set_ydata([a_["thr_db"]] * 2)

    good = (a_["snr"] >= MIN_SNR_DB and a_["n_on"] >= MIN_SAMPLES and a_["coh"] >= MIN_COH)
    if good:
        dphi_raw = np.angle(a_["corr"])
        if state["do_cal"]:
            state["cal"], state["do_cal"] = dphi_raw, False
            trail.clear(); hist.clear(); est = None
            print(f"[cal] offset {np.degrees(dphi_raw):+.1f} deg")
        dphi = np.angle(np.exp(1j * (dphi_raw - state["cal"])))
        sgn = state["sign"]
        hist.append(np.arcsin(np.clip(-sgn * dphi / k_d, -1, 1)))
        med = float(np.median(hist))                         # outliers never reach the needle
        est = med if est is None else ALPHA * est + (1 - ALPHA) * med

        S = a_["S"].copy(); S[1] *= np.exp(-1j * state["cal"])
        a = np.vstack([np.ones_like(th), np.exp(-1j * sgn * k_d * np.sin(th))])
        R = S @ S.conj().T / S.shape[1]
        if state["music"]:
            _, V = np.linalg.eigh(R)
            P = 1.0 / (np.abs(V[:, :1].conj().T @ a).ravel() ** 2 + 1e-12)
        else:
            P = np.real(np.sum(a.conj() * (R @ a), axis=0))
        r = (np.clip(10 * np.log10(P / P.max()), -30, 0) + 30) / 30
        beam_line.set_data(th, r)
        beam_fill.set_xy(np.column_stack([np.r_[th, th[::-1]], np.r_[r, np.zeros_like(r)]]))

        needle.set_data([est, est], [0, 1]); glow.set_data([est, est], [0, 1])
        needle.set_alpha(1.0)
        trail.append(est)
        n = len(trail)
        trail_sc.set_offsets(np.column_stack([list(trail), np.linspace(0.55, 0.95, n)]))
        trail_sc.set_array(np.linspace(0, 1, n))
        trail_sc.set_sizes(np.linspace(5, 60, n))
        big.set_text(f"{np.degrees(est):+6.1f}°"); big.set_color(MAG)
        ph1.set_data([0, np.cos(dphi)], [0, np.sin(dphi)])
        ph_txt.set_text(f"Δφ = {np.degrees(dphi):+6.1f}°")
    else:
        if state["do_cal"]:
            print("[cal] no clean signal - calibration skipped")
            state["do_cal"] = False
        needle.set_alpha(0.3)                                  # hold last bearing, dimmed
        if a_["snr"] < MIN_SNR_DB:
            big.set_text("NO SIGNAL"); big.set_color(RED)

    now = time.time()
    fps = 0.9 * fps + 0.1 / max(now - t_prev, 1e-3); t_prev = now
    clip = a_["pk"] > 0.95
    info.set_text(f"{'MUSIC' if state['music'] else 'Bartlett'} | {a_['mode']} | SNR {a_['snr']:4.1f} dB | "
                  f"coh {a_['coh']:.2f} | gain {state['gain']} dB | {fps:4.1f} fps"
                  + ("  | CLIPPING - press g" if clip else ""))
    info.set_color(RED if clip else GREEN)

    if now - t_log > 1.0:
        t_log = now
        print(f"[diag] NB {a_['snr_nb']:5.1f} dB  TD {a_['snr_td']:5.1f} dB  used {a_['mode']}  "
              f"coh {a_['coh']:.2f}  samples {a_['n_on']:6d}  peak {a_['pk']*100:5.1f}% FS  "
              f"line {a_['f_pk']:+7.1f} kHz  gain {state['gain']}  {fps:4.1f} fps")
    render()
