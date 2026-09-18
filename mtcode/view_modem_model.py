#!/usr/bin/env python3
"""
Interactive top-view viewer and editor for ModEM 3-D resistivity models.

Edit the settings at the top, then run this file (Run / F5).

1. Click 4 points on the map
2. Set layer range and resistivity in the boxes below the map
3. Click the green "Apply + Save" button (native Windows button, not inside the plot)
"""

from __future__ import annotations

import sys
import traceback
from pathlib import Path
from tkinter import (
    BOTH,
    BOTTOM,
    LEFT,
    TOP,
    X,
    Button,
    Entry,
    Frame,
    Label,
    Scale,
    StringVar,
    Tk,
    messagebox,
)

import matplotlib

matplotlib.use("TkAgg")

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import Polygon

from modem_model import ModemModel, load_modem_model, save_modem_model

# --- user settings: edit these, then click Run ---
SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_PATH = SCRIPT_DIR.parent / "Forward/Dike1/2mgrid/Modelfile.model"
OUTPUT_PATH = SCRIPT_DIR.parent / "Forward/Dike1/2mgrid/MyEditedModel.model"
START_LAYER = 1
DEFAULT_LAYER_FROM = 5
DEFAULT_LAYER_TO = 15
DEFAULT_RESISTIVITY = 28.0
NUM_PICK_POINTS = 4
# -------------------------------------------------


def _output_path_for(model_path: Path) -> Path:
    if OUTPUT_PATH is not None:
        return Path(OUTPUT_PATH).resolve()
    return model_path.resolve().with_name(f"{model_path.stem}_edited{model_path.suffix}")


def _jet_like_cmap():
    base = plt.cm.jet(np.linspace(0, 1, 24))
    keep = [0, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
    return mcolors.ListedColormap(base[keep][::-1])


def plot_layer(ax, model, layer, vmin, vmax, mesh, title):
    data = model.layer_slice(layer)
    log_data = np.log10(data)

    y_km = model.y_edges / 1000.0
    x_km = model.x_edges / 1000.0
    yy, xx = np.meshgrid(y_km, x_km)

    if mesh is None:
        mesh = ax.pcolormesh(
            yy,
            xx,
            log_data,
            shading="flat",
            cmap=_jet_like_cmap(),
            norm=plt.Normalize(vmin=np.log10(vmin), vmax=np.log10(vmax)),
        )
    else:
        mesh.set_array(log_data.ravel())
        mesh.set_clim(np.log10(vmin), np.log10(vmax))

    z0, z1 = model.layer_depth_range_km(layer)
    title.set_text(
        f"Layer {layer + 1}/{model.nz} | depth {z0:.3f} to {z1:.3f} km"
    )
    ax.set_xlabel("Easting Y (km)")
    ax.set_ylabel("Northing X (km)")
    ax.set_aspect("equal", adjustable="box")
    return mesh


def _finite_resistivity_bounds(model: ModemModel) -> tuple[float, float]:
    finite = model.rho[np.isfinite(model.rho) & (model.rho > 0)]
    if finite.size == 0:
        return 1.0, 1000.0
    return float(np.min(finite)), float(np.max(finite))


def interactive_view(
    model: ModemModel,
    start_layer: int = 1,
    output_path: Path | None = None,
) -> None:
    start_layer = int(np.clip(start_layer, 1, model.nz))
    layer_idx = start_layer - 1
    vmin, vmax = _finite_resistivity_bounds(model)

    root = Tk()
    root.title("ModEM model editor")
    root.minsize(900, 700)
    root.geometry("1050x850")

    # Controls pinned to bottom so the plot cannot push them off-screen.
    panel = Frame(root, padx=10, pady=8, relief="ridge", bd=1)
    panel.pack(side=BOTTOM, fill=X)

    fig = Figure(figsize=(9, 6), dpi=100)
    ax = fig.add_subplot(111)
    title = ax.set_title("")
    mesh = plot_layer(ax, model, layer_idx, vmin, vmax, None, title)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    tick_vals = [0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
    tick_vals = [v for v in tick_vals if vmin <= v <= vmax] or [vmin, vmax]
    cbar.set_ticks(np.log10(tick_vals))
    cbar.set_ticklabels([str(v) for v in tick_vals])
    cbar.set_label("Resistivity (ohm-m)")
    fig.suptitle(f"{model.nx}x{model.ny}x{model.nz} cells  |  saves to: {output_path.name}")

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True, padx=8, pady=8)

    status_var = StringVar(
        value=f"Step 1: click {NUM_PICK_POINTS} points on the map."
    )
    status_label = Label(
        panel,
        textvariable=status_var,
        anchor="w",
        justify="left",
        bg="#eef3ff",
        fg="#111",
        font=("Segoe UI", 10),
        padx=8,
        pady=6,
        relief="groove",
    )
    status_label.pack(fill=X, pady=(0, 6))

    slider_row = Frame(panel)
    slider_row.pack(fill=X, pady=4)

    layer_label_var = StringVar(value=f"View layer: {start_layer} / {model.nz}")
    Label(slider_row, textvariable=layer_label_var, font=("Segoe UI", 10), width=18, anchor="w").pack(
        side=LEFT, padx=(0, 8)
    )
    layer_scale = Scale(
        slider_row,
        from_=1,
        to=model.nz,
        orient="horizontal",
        length=500,
        resolution=1,
        showvalue=True,
    )
    layer_scale.set(start_layer)
    layer_scale.pack(side=LEFT, fill=X, expand=True)

    edit_row = Frame(panel)
    edit_row.pack(fill=X, pady=4)

    Label(edit_row, text="Edit layers — from:", font=("Segoe UI", 10)).pack(side=LEFT)
    entry_from = Entry(edit_row, width=6, font=("Segoe UI", 11), relief="solid", bd=1)
    entry_from.insert(0, str(DEFAULT_LAYER_FROM))
    entry_from.pack(side=LEFT, padx=(6, 16))

    Label(edit_row, text="to:", font=("Segoe UI", 10)).pack(side=LEFT)
    entry_to = Entry(edit_row, width=6, font=("Segoe UI", 11), relief="solid", bd=1)
    entry_to.insert(0, str(DEFAULT_LAYER_TO))
    entry_to.pack(side=LEFT, padx=(6, 16))

    Label(edit_row, text="Rho (ohm-m):", font=("Segoe UI", 10)).pack(side=LEFT)
    entry_rho = Entry(edit_row, width=10, font=("Segoe UI", 11), relief="solid", bd=1)
    entry_rho.insert(0, str(DEFAULT_RESISTIVITY))
    entry_rho.pack(side=LEFT, padx=(6, 0))

    state = {
        "points": [],
        "markers": [],
        "polygon": None,
        "dirty": False,
    }

    def set_status(msg: str) -> None:
        status_var.set(msg)
        print(msg)
        root.update_idletasks()

    def redraw_current_layer() -> None:
        plot_layer(ax, model, layer_idx, vmin, vmax, mesh, title)
        canvas.draw_idle()

    def on_layer_scale(_value) -> None:
        nonlocal layer_idx
        layer_idx = int(float(layer_scale.get())) - 1
        layer_label_var.set(f"View layer: {layer_idx + 1} / {model.nz}")
        redraw_current_layer()

    layer_scale.configure(command=on_layer_scale)

    def clear_picks() -> None:
        state["points"] = []
        for marker in state["markers"]:
            marker.remove()
        state["markers"] = []
        if state["polygon"] is not None:
            state["polygon"].remove()
            state["polygon"] = None
        canvas.draw_idle()

    def update_polygon_outline() -> None:
        if state["polygon"] is not None:
            state["polygon"].remove()
            state["polygon"] = None
        if len(state["points"]) < 2:
            return
        pts = state["points"]
        closed = pts + [pts[0]] if len(pts) == NUM_PICK_POINTS else pts
        state["polygon"] = ax.add_patch(
            Polygon(
                closed,
                closed=False,
                linewidth=2,
                edgecolor="yellow",
                facecolor="yellow",
                alpha=0.2,
            )
        )

    def on_map_click(event) -> None:
        if event.inaxes is not ax or event.button != 1:
            return
        if event.xdata is None or event.ydata is None:
            return

        if len(state["points"]) >= NUM_PICK_POINTS:
            clear_picks()

        y_km, x_km = float(event.xdata), float(event.ydata)
        state["points"].append((y_km, x_km))
        marker = ax.plot(
            y_km, x_km, "wo", markeredgecolor="black", markersize=8, zorder=5
        )[0]
        state["markers"].append(marker)
        update_polygon_outline()
        canvas.draw_idle()

        n = len(state["points"])
        if n < NUM_PICK_POINTS:
            set_status(f"Step 1: point {n}/{NUM_PICK_POINTS} placed. Click next corner.")
        else:
            set_status(
                f"Step 2: {NUM_PICK_POINTS} points set. "
                "Edit layers/rho, then click the green Apply + Save button below."
            )

    canvas.mpl_connect("button_press_event", on_map_click)

    def parse_int(text: str, name: str) -> int:
        value = int(float(text.strip()))
        if value < 1 or value > model.nz:
            raise ValueError(f"{name} must be between 1 and {model.nz}")
        return value

    def on_apply() -> None:
        set_status("Working: apply + save ...")
        try:
            if len(state["points"]) < NUM_PICK_POINTS:
                messagebox.showwarning(
                    "Need 4 points",
                    f"Click {NUM_PICK_POINTS} points on the map first.\n"
                    f"You have placed {len(state['points'])}.",
                )
                set_status(f"Need {NUM_PICK_POINTS} points on the map (have {len(state['points'])}).")
                return

            layer_from = parse_int(entry_from.get(), "Layer from")
            layer_to = parse_int(entry_to.get(), "Layer to")
            resistivity = float(entry_rho.get().strip())
            if resistivity <= 0:
                raise ValueError("Resistivity must be positive")

            vertices = np.array(state["points"], dtype=float)
            in_polygon = int(model.mask_cells_in_polygon_km(vertices).sum())

            changed = model.apply_resistivity_polygon(
                vertices, layer_from, layer_to, resistivity
            )
            if changed == 0:
                messagebox.showwarning(
                    "No cells updated",
                    f"{in_polygon} cells inside your area, but layers "
                    f"{layer_from}-{layer_to} are air there.\n\n"
                    "Try deeper layers (e.g. from 10 to 20).",
                )
                set_status("No earth cells in that layer range. Try deeper layers.")
                return

            state["dirty"] = True
            redraw_current_layer()

            set_status(f"Saving to {output_path} ...")
            root.update_idletasks()
            saved = save_modem_model(model, output_path)
            state["dirty"] = False

            msg = (
                f"Updated {changed} cells.\n"
                f"Rho = {resistivity} ohm-m, layers {layer_from}-{layer_to}.\n\n"
                f"Saved to:\n{saved}"
            )
            set_status(f"Done. Saved {changed} cells to {saved.name}")
            messagebox.showinfo("Success", msg)
        except Exception as exc:
            set_status(f"Error: {exc}")
            messagebox.showerror("Failed", f"{exc}\n\n{traceback.format_exc()}")

    def on_clear() -> None:
        clear_picks()
        set_status(f"Points cleared. Click {NUM_PICK_POINTS} points on the map.")

    btn_row = Frame(panel)
    btn_row.pack(fill=X, pady=(8, 0))

    apply_btn = Button(
        btn_row,
        text="Apply + Save",
        command=on_apply,
        font=("Segoe UI", 11, "bold"),
        bg="#2e9b3c",
        fg="white",
        activebackground="#248a31",
        activeforeground="white",
        padx=16,
        pady=8,
        cursor="hand2",
    )
    apply_btn.pack(side=LEFT, padx=(0, 8))

    Button(
        btn_row,
        text="Clear points",
        command=on_clear,
        font=("Segoe UI", 10),
        padx=12,
        pady=8,
    ).pack(side=LEFT)

    def on_close() -> None:
        if state["dirty"]:
            try:
                save_modem_model(model, output_path)
                print(f"Saved on close: {output_path}")
            except Exception as exc:
                print(f"Save on close failed: {exc}", file=sys.stderr)
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)
    root.bind("<Return>", lambda _e: on_apply())
    root.mainloop()


def main() -> int:
    model_path = Path(MODEL_PATH)
    if not model_path.is_file():
        print(f"Model file not found: {model_path}", file=sys.stderr)
        return 1

    output_path = _output_path_for(model_path)
    model = load_modem_model(model_path)
    print(f"Loaded: {model_path.resolve()}")
    print(f"Output: {output_path}")
    interactive_view(model, start_layer=START_LAYER, output_path=output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
