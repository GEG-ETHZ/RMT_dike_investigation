"""Read and write ModEM 3-D resistivity models (WSINV / M3D format)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


AIR_RESISTIVITY = 1e17


@dataclass
class ModemModel:
    """In-memory ModEM model aligned with mtcode ``load_model_modem.m``."""

    nx: int
    ny: int
    nz: int
    nz_air: int
    rho_type: str
    dx: np.ndarray
    dy: np.ndarray
    dz: np.ndarray
    rho: np.ndarray  # shape (nx, ny, nz), linear resistivity in ohm-m
    origin: np.ndarray  # [south, west, top] in metres
    rotation: int
    header_lines: list[str]
    tail_lines: list[str]

    @property
    def x_edges(self) -> np.ndarray:
        return np.concatenate(([0.0], np.cumsum(self.dx))) + self.origin[0]

    @property
    def y_edges(self) -> np.ndarray:
        return np.concatenate(([0.0], np.cumsum(self.dy))) + self.origin[1]

    @property
    def z_edges(self) -> np.ndarray:
        return np.concatenate(([0.0], np.cumsum(self.dz))) + self.origin[2]

    @property
    def cx(self) -> np.ndarray:
        x = self.x_edges
        return (x[:-1] + x[1:]) / 2.0

    @property
    def cy(self) -> np.ndarray:
        y = self.y_edges
        return (y[:-1] + y[1:]) / 2.0

    @property
    def cz(self) -> np.ndarray:
        z = self.z_edges
        return (z[:-1] + z[1:]) / 2.0

    def layer_slice(self, layer: int) -> np.ndarray:
        """Return resistivity map for one z-layer (nx, ny), NaN for air cells."""
        data = self.rho[:, :, layer].astype(float).copy()
        data[data > 1e15] = np.nan
        return data

    def layer_depth_range_km(self, layer: int) -> tuple[float, float]:
        z = self.z_edges / 1000.0
        return float(z[layer]), float(z[layer + 1])

    def mask_cells_in_box_km(
        self,
        y_km_min: float,
        y_km_max: float,
        x_km_min: float,
        x_km_max: float,
    ) -> np.ndarray:
        """
        Boolean mask (nx, ny) for cell centres inside a plot-coordinate box.

        Plot axes match M3D: horizontal = easting Y, vertical = northing X.
        """
        y0, y1 = sorted((y_km_min, y_km_max))
        x0, x1 = sorted((x_km_min, x_km_max))
        x_in = (self.cx / 1000.0 >= x0) & (self.cx / 1000.0 <= x1)
        y_in = (self.cy / 1000.0 >= y0) & (self.cy / 1000.0 <= y1)
        return x_in[:, None] & y_in[None, :]

    def mask_cells_in_polygon_km(self, vertices_yx_km: np.ndarray) -> np.ndarray:
        """
        Boolean mask (nx, ny) for cell centres inside a polygon.

        ``vertices_yx_km`` is (n, 2) with columns [easting Y km, northing X km].
        """
        from matplotlib.path import Path

        if vertices_yx_km.shape[0] < 3:
            return np.zeros((self.nx, self.ny), dtype=bool)

        y_plot, x_plot = np.meshgrid(self.cy / 1000.0, self.cx / 1000.0)
        points = np.column_stack([y_plot.ravel(), x_plot.ravel()])
        inside = Path(vertices_yx_km).contains_points(points)
        return inside.reshape(self.nx, self.ny)

    def _apply_resistivity_mask_xy(
        self,
        mask_xy: np.ndarray,
        layer_from: int,
        layer_to: int,
        resistivity: float,
        *,
        skip_air: bool = True,
    ) -> int:
        if layer_from > layer_to:
            layer_from, layer_to = layer_to, layer_from
        layer_from = int(np.clip(layer_from, 1, self.nz))
        layer_to = int(np.clip(layer_to, 1, self.nz))

        if not mask_xy.any():
            return 0

        ii, jj = np.where(mask_xy)
        changed = 0
        for k in range(layer_from - 1, layer_to):
            slab = self.rho[ii, jj, k]
            if skip_air:
                earth = np.isfinite(slab)
                self.rho[ii[earth], jj[earth], k] = resistivity
                changed += int(earth.sum())
            else:
                self.rho[ii, jj, k] = resistivity
                changed += len(ii)
        return changed

    def apply_resistivity_polygon(
        self,
        vertices_yx_km: np.ndarray,
        layer_from: int,
        layer_to: int,
        resistivity: float,
        *,
        skip_air: bool = True,
    ) -> int:
        """Set resistivity inside a polygon; layers are 1-based inclusive."""
        mask_xy = self.mask_cells_in_polygon_km(vertices_yx_km)
        return self._apply_resistivity_mask_xy(
            mask_xy, layer_from, layer_to, resistivity, skip_air=skip_air
        )

    def apply_resistivity_box(
        self,
        y_km_min: float,
        y_km_max: float,
        x_km_min: float,
        x_km_max: float,
        layer_from: int,
        layer_to: int,
        resistivity: float,
        *,
        skip_air: bool = True,
    ) -> int:
        """
        Set resistivity inside a horizontal box over a layer range.

        ``layer_from`` and ``layer_to`` are 1-based and inclusive.
        Returns the number of cells updated.
        """
        y0, y1 = sorted((y_km_min, y_km_max))
        x0, x1 = sorted((x_km_min, x_km_max))
        mask_xy = self.mask_cells_in_box_km(y0, y1, x0, x1)
        return self._apply_resistivity_mask_xy(
            mask_xy, layer_from, layer_to, resistivity, skip_air=skip_air
        )


def load_modem_model(path: str | Path) -> ModemModel:
    """
    Load a ModEM model file.

    Indexing matches ``load_model_modem.m`` / ``read_mackie3d_model`` with
    ``block=true``: ``rho[i, j, k]`` is north-south index ``i``, east-west ``j``,
    depth layer ``k`` (0-based).
    """
    lines = Path(path).read_text().splitlines()
    header: list[str] = []
    idx = 0
    while idx < len(lines) and (
        not lines[idx].strip()
        or lines[idx].lstrip().startswith("#")
    ):
        header.append(lines[idx])
        idx += 1

    dims = lines[idx].split()
    nx, ny, nz, nz_air = map(int, dims[:4])
    rho_type = dims[4] if len(dims) > 4 else "LINEAR"
    idx += 1

    dx = np.array([float(v) for v in lines[idx].split()], dtype=float)
    idx += 1
    dy = np.array([float(v) for v in lines[idx].split()], dtype=float)
    idx += 1
    dz = np.array([float(v) for v in lines[idx].split()], dtype=float)
    idx += 1

    if len(dx) != nx or len(dy) != ny or len(dz) != nz:
        raise ValueError(
            f"Grid size mismatch in {path}: header ({nx}, {ny}, {nz}) vs "
            f"spacing vectors ({len(dx)}, {len(dy)}, {len(dz)})"
        )

    # Skip blank line before first layer, if present.
    while idx < len(lines) and not lines[idx].strip():
        idx += 1

    raw_layers: list[np.ndarray] = []
    for _ in range(nz):
        rows = []
        for _ in range(ny):
            while idx < len(lines) and not lines[idx].strip():
                idx += 1
            if idx >= len(lines):
                raise ValueError(f"Unexpected end of file while reading layer data in {path}")
            rows.append([float(v) for v in lines[idx].split()])
            idx += 1
        layer = np.array(rows, dtype=float)  # (ny, nx), file row = y, reversed x
        raw_layers.append(layer)

    tail = lines[idx:]
    origin = np.array([-sum(dx) / 2.0, -sum(dy) / 2.0, 0.0], dtype=float)
    rotation = 0
    for line in tail:
        values = [float(v) for v in line.split()]
        if len(values) == 3:
            origin = np.array(values, dtype=float)
        elif len(values) == 1:
            rotation = int(values[0])

    # Match MATLAB ``flipud`` applied to the (nx, ny) block read from file.
    rho = np.stack([np.flipud(layer.T) for layer in raw_layers], axis=2)

    if rho_type.upper() == "LOGE":
        rho = np.exp(rho)

    rho[rho > 1e15] = np.nan

    return ModemModel(
        nx=nx,
        ny=ny,
        nz=nz,
        nz_air=nz_air,
        rho_type=rho_type,
        dx=dx,
        dy=dy,
        dz=dz,
        rho=rho,
        origin=origin,
        rotation=rotation,
        header_lines=header,
        tail_lines=tail,
    )


def save_modem_model(model: ModemModel, path: str | Path) -> Path:
    """Write model back to ModEM format. Returns the resolved output path."""
    path = Path(path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)

    rho_out = model.rho.copy()
    rho_out[np.isnan(rho_out)] = AIR_RESISTIVITY

    if model.rho_type.upper() == "LOGE":
        rho_out = np.log(rho_out)

    with path.open("w", encoding="utf-8", newline="\n") as fid:
        if model.header_lines:
            for line in model.header_lines:
                fid.write(f"{line}\n")
        else:
            fid.write("# Written by Python modem_model.py\n")

        fid.write(
            f"{model.nx} {model.ny} {model.nz} {model.nz_air} {model.rho_type}\n"
        )
        fid.write(" ".join(f"{v:.6f}" for v in model.dx) + "\n")
        fid.write(" ".join(f"{v:.6f}" for v in model.dy) + "\n")
        fid.write(" ".join(f"{v:.6f}" for v in model.dz) + "\n")

        for k in range(model.nz):
            fid.write("\n")
            for j in range(model.ny):
                # Match MATLAB write_model_modem: i = nx:-1:1 on each row.
                row = rho_out[::-1, j, k]
                fid.write("    " + "    ".join(f"{v:.5E}" for v in row) + "\n")

        if model.tail_lines:
            for line in model.tail_lines:
                fid.write(f"{line}\n")
        else:
            fid.write(
                f"{model.origin[0]:.4f} {model.origin[1]:.4f} "
                f"{model.origin[2]:.4f}\n"
            )
            fid.write(f"{model.rotation}\n")

    if not path.is_file() or path.stat().st_size == 0:
        raise OSError(f"Save failed: file was not written: {path}")

    return path
