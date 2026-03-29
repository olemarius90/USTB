"""Plotting utilities for USTB beamformed data."""

import numpy as np
import matplotlib.pyplot as plt


def plot_beamformed_data(b_data, title="", dynamic_range=60):
    """Plot beamformed data as a B-mode image.

    Mirrors the MATLAB beamformed_data.plot() method.
    """
    scan = b_data.scan
    data = np.array(b_data.data)

    if data.ndim == 4:
        data = data[:, 0, 0, 0]
    elif data.ndim == 3:
        data = data[:, 0, 0]
    elif data.ndim == 2:
        data = data[:, 0]

    envelope = np.abs(data)
    envelope[envelope < np.finfo(float).eps] = np.finfo(float).eps
    img_db = 20.0 * np.log10(envelope / envelope.max())
    img_db = np.clip(img_db, -dynamic_range, 0)

    if hasattr(scan, "azimuth_axis") and scan.azimuth_axis is not None:
        N_az = len(scan.azimuth_axis)
        N_depth = len(scan.depth_axis)
        img_2d = img_db.reshape(N_az, N_depth).T

        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        extent = [
            np.rad2deg(scan.azimuth_axis[0]),
            np.rad2deg(scan.azimuth_axis[-1]),
            scan.depth_axis[-1] * 1e3,
            scan.depth_axis[0] * 1e3,
        ]
        ax.imshow(img_2d, cmap="gray", aspect="auto", extent=extent,
                  vmin=-dynamic_range, vmax=0)
        ax.set_xlabel("Azimuth [deg]")
        ax.set_ylabel("Depth [mm]")
    elif hasattr(scan, "x_axis") and scan.x_axis is not None:
        N_x = len(scan.x_axis)
        N_z = len(scan.z_axis)
        img_2d = img_db.reshape(N_x, N_z).T

        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        extent = [
            scan.x_axis[0] * 1e3,
            scan.x_axis[-1] * 1e3,
            scan.z_axis[-1] * 1e3,
            scan.z_axis[0] * 1e3,
        ]
        ax.imshow(img_2d, cmap="gray", aspect="auto", extent=extent,
                  vmin=-dynamic_range, vmax=0)
        ax.set_xlabel("Lateral [mm]")
        ax.set_ylabel("Depth [mm]")
    else:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        ax.plot(img_db)
        ax.set_ylabel("Amplitude [dB]")

    ax.set_title(title if title else "B-mode Image")
    plt.colorbar(ax.images[0] if ax.images else None, ax=ax, label="dB")
    plt.tight_layout()
    return fig, ax
