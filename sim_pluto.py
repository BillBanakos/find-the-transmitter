"""
Simulated PlutoSDR+ for `python pluto_aoa_radar.py --sim` (no hardware needed).

Mimics the parts of pyadi-iio's `adi.ad9361` the demo uses. A virtual transmitter sends
short noise-like bursts (like Wi-Fi beacons or SA pulses) that reach two antennas spaced
SPACING metres apart. It stands at 0 deg for the first CAL_TIME seconds (press c then),
then walks back and forth across the arc. A weaker interferer sits at a fixed angle.
"""
import time
import numpy as np

SPACING = 0.0625          # m, same as D in pluto_aoa_radar.py
CAL_TIME = 6.0            # s at 0 deg before walking starts
WALK_DEG, WALK_PERIOD = 50.0, 12.0
BURST, PRI = 0.4e-3, 1.0e-3
CH_OFFSET = 0.7           # rad, random RX1/RX2 phase offset of a real AD9361 (removed by 'c')


class _RxAdc:
    def set_kernel_buffers_count(self, n):
        pass


class ad9361:
    def __init__(self, uri="sim"):
        self._rxadc = _RxAdc()
        self.rx_hardwaregain_chan0 = 40
        self.rx_hardwaregain_chan1 = 40
        self.rx_lo = int(2.478e9)
        self.sample_rate = int(2e6)
        self.rx_buffer_size = 2**15
        self._t0 = time.time()
        self._rng = np.random.default_rng()

    def true_angle(self):
        t = time.time() - self._t0
        if t < CAL_TIME:
            return 0.0
        return WALK_DEG * np.sin(2 * np.pi * (t - CAL_TIME) / WALK_PERIOD)

    def rx(self):
        n, fs, rng = self.rx_buffer_size, self.sample_rate, self._rng
        kd = 2 * np.pi * SPACING / (3e8 / self.rx_lo)
        g = 10 ** ((self.rx_hardwaregain_chan0 - 40) / 20)
        noise = lambda: 3 * g * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
        x0, x1 = noise(), noise()
        L = int(BURST * fs)

        def add_bursts(angle_deg, amp, pri):
            start = rng.uniform(0, pri)
            while int(start * fs) + L < n:
                i = int(start * fs)
                s = amp * g * (rng.standard_normal(L) + 1j * rng.standard_normal(L)) / np.sqrt(2)
                ph = kd * np.sin(np.radians(angle_deg)) + np.radians(3) * rng.standard_normal()
                x0[i:i + L] += s
                x1[i:i + L] += s * np.exp(-1j * ph) * np.exp(1j * CH_OFFSET)
                start += pri

        add_bursts(self.true_angle(), 300, PRI)                 # our transmitter
        add_bursts(-60, 300 * 10 ** (-15 / 20), PRI * 7.3)      # weaker interferer
        clip = lambda x: np.clip(x.real, -2048, 2047) + 1j * np.clip(x.imag, -2048, 2047)
        time.sleep(n / fs)                                      # real-time pacing
        return [clip(x0 + 40), clip(x1 + 40)]                   # + DC (LO leakage)
