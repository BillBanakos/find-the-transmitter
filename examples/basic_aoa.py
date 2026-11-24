"""
"Find the transmitter" - live 2-element AoA demo on PlutoSDR+ (AD9361, 2R2T).

Setup
  RX1, RX2 : two identical monopoles, spacing D = lambda/2, both vertical.
  TX1      : third monopole on a long SMA cable, held by a volunteer who walks around the array.
  Pluto+ must be in 2r2t mode (once, over ssh):
      fw_setenv attr_name compatible; fw_setenv attr_val ad9361; fw_setenv mode 2r2t; reboot

Run:  python pluto_aoa_demo.py
Keys: c = calibrate (volunteer standing at broadside, 0 deg)   f = flip direction   q = quit
Recalibrate after any LO retune or reboot: the RX1/RX2 phase offset is not fixed across retunes.
"""
import numpy as np
import adi
import matplotlib.pyplot as plt

URI = "ip:192.168.2.1"
F_SIG = 2.478e9        # SA tone, inside the 2400-2483.5 MHz SRD band
FC = F_SIG - 200e3     # Pluto LO 200 kHz below -> tone lands at +200 kHz, away from DC leakage
FS = 2.0e6
N = 2**14
LAM = 3e8 / F_SIG      # 0.125 m -> lambda/2 = 6.25 cm
D = LAM / 2            # element spacing [m]; measure the real value and put it here
TONE = FS * 1600 / N   # integer number of cycles per buffer -> clean cyclic TX
GAIN = 40              # manual RX gain [dB]; lower it if the tone clips

sdr = adi.ad9361(uri=URI)
sdr.sample_rate = int(FS)
sdr.rx_lo = int(FC)
sdr.tx_lo = int(FC)
sdr.rx_rf_bandwidth = int(FS)
sdr.tx_rf_bandwidth = int(FS)
sdr.rx_enabled_channels = [0, 1]
sdr.rx_buffer_size = N
sdr.gain_control_mode_chan0 = "manual"
sdr.gain_control_mode_chan1 = "manual"
sdr.rx_hardwaregain_chan0 = GAIN
sdr.rx_hardwaregain_chan1 = GAIN
USE_PLUTO_TX = False   # False: external CW source (e.g. SA tracking generator, zero span, at FC + ~200 kHz)
if USE_PLUTO_TX:
    sdr.tx_enabled_channels = [0]
    sdr.tx_hardwaregain_chan0 = -10
    sdr.tx_cyclic_buffer = True
    t = np.arange(N) / FS
    sdr.tx(0.5 * 2**14 * np.exp(2j * np.pi * TONE * t))
else:
    sdr.tx_hardwaregain_chan0 = -89                # keep Pluto's own TX/LO leakage quiet
    sdr.tx_hardwaregain_chan1 = -89

win = np.hanning(N)
DC_GUARD = int(20e3 / FS * N)                      # ignore +-20 kHz around DC (LO leakage)
k_d = 2 * np.pi * D / LAM
thetas = np.radians(np.linspace(-90, 90, 361))
state = {"cal": 0.0, "sign": 1.0, "run": True}


def snapshot():
    x0, x1 = sdr.rx()
    X0 = np.fft.fft(x0 * win)
    X1 = np.fft.fft(x1 * win)
    mag = np.abs(X0) + np.abs(X1)
    mag[:DC_GUARD] = 0; mag[-DC_GUARD:] = 0
    k = np.argmax(mag)                               # strongest bin = our tone
    v = np.array([X0[k], X1[k]])
    snr_db = 20 * np.log10(np.abs(X0[k]) / np.median(np.abs(X0)))
    dphi = np.angle(v[1] * np.conj(v[0]))
    return v, dphi, snr_db


def on_key(ev):
    if ev.key == "c":
        _, state["cal"], _ = snapshot()
        print(f"calibrated: offset = {np.degrees(state['cal']):.1f} deg")
    elif ev.key == "f":
        state["sign"] *= -1
    elif ev.key == "q":
        state["run"] = False


fig = plt.figure(figsize=(7, 5))
ax = fig.add_subplot(projection="polar")
ax.set_thetamin(-90); ax.set_thetamax(90); ax.set_theta_zero_location("N")
spec_line, = ax.plot(thetas, np.zeros_like(thetas), lw=2)
needle, = ax.plot([0, 0], [0, 1], "r", lw=3)
title = ax.set_title("")
fig.canvas.mpl_connect("key_press_event", on_key)

ema = None
while state["run"] and plt.fignum_exists(fig.number):
    v, dphi, snr = snapshot()
    v[1] *= np.exp(-1j * state["cal"])                # remove channel phase offset
    dphi = np.angle(np.exp(1j * (dphi - state["cal"])))
    s = np.clip(-state["sign"] * dphi / k_d, -1, 1)
    theta = np.arcsin(s)
    ema = theta if ema is None else 0.7 * ema + 0.3 * theta

    a = np.vstack([np.ones_like(thetas), np.exp(-1j * state["sign"] * k_d * np.sin(thetas))])
    P = np.abs(a.conj().T @ v) ** 2                  # Bartlett (delay-and-sum) scan
    spec_line.set_ydata(P / P.max())
    needle.set_xdata([ema, ema])
    title.set_text(f"AoA = {np.degrees(ema):+5.1f} deg   dphi = {np.degrees(dphi):+6.1f} deg   SNR ~ {snr:4.1f} dB")
    plt.pause(0.05)

if USE_PLUTO_TX:
    sdr.tx_destroy_buffer()
