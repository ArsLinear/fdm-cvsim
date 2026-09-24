import numpy as np

def cv_waveform_generator(initial, switch, scan_rate, dt):

    t_1 = np.arange(0, (switch - initial) / scan_rate, dt)
    t_2 = np.arange((switch - initial) / scan_rate, 2 * (switch - initial) / scan_rate, dt)

    waveform_1 = initial + scan_rate * t_1
    waveform_2 = switch - scan_rate * (t_2 - (switch - initial) / scan_rate)

    waveform = np.concatenate([waveform_1, waveform_2])

    return waveform

def ftacv_waveform_generator(initial, switch, scan_rate, amplitude, frequency, dt):

    t_1 = np.arange(0, (switch - initial) / scan_rate, dt)
    t_2 = np.arange((switch - initial) / scan_rate, 2 * (switch - initial) / scan_rate, dt)

    waveform_1 = initial + scan_rate * t_1 + amplitude * np.sin(2 * np.pi * frequency * t_1)
    waveform_2 = switch - scan_rate * (t_2 - (switch - initial) / scan_rate) + amplitude * np.sin(2 * np.pi * frequency * t_2)

    waveform = np.concatenate([waveform_1, waveform_2])

    return waveform