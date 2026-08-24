"""
ui/widgets.py
=============

Small reusable tkinter widgets:

* :class:`ScrollableFrame` — a vertically scrollable container (Canvas +
  Scrollbar) so a tall panel of controls never gets clipped on any display.
  Put content into ``.inner``; grid/pack ``.outer``.
* :class:`Tooltip` — a lightweight hover tooltip for any widget.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk


class ScrollableFrame:
    """A vertically scrollable frame.

    Usage::

        sf = ScrollableFrame(parent)
        sf.outer.grid(row=0, column=1, sticky="nsew")
        something = ttk.Label(sf.inner, text="...")   # add content to .inner
    """

    def __init__(self, parent, bg: str = "#252526") -> None:
        self.outer = ttk.Frame(parent)
        self.outer.rowconfigure(0, weight=1)
        self.outer.columnconfigure(0, weight=1)

        self._canvas = tk.Canvas(self.outer, borderwidth=0, highlightthickness=0,
                                 background=bg)
        self._vsb = ttk.Scrollbar(self.outer, orient="vertical",
                                  command=self._canvas.yview)
        self._canvas.configure(yscrollcommand=self._vsb.set)
        self._canvas.grid(row=0, column=0, sticky="nsew")
        self._vsb.grid(row=0, column=1, sticky="ns")

        self.inner = ttk.Frame(self._canvas)
        self._win = self._canvas.create_window((0, 0), window=self.inner,
                                               anchor="nw")

        self.inner.bind("<Configure>", self._on_inner_configure)
        self._canvas.bind("<Configure>", self._on_canvas_configure)
        # Wheel scrolling only while the pointer is over this canvas.
        self._canvas.bind("<Enter>", self._bind_wheel)
        self._canvas.bind("<Leave>", self._unbind_wheel)

    def _on_inner_configure(self, _event) -> None:
        self._canvas.configure(scrollregion=self._canvas.bbox("all"))

    def _on_canvas_configure(self, event) -> None:
        # Make the inner frame track the canvas width (so widgets fill across).
        self._canvas.itemconfigure(self._win, width=event.width)

    def _bind_wheel(self, _event) -> None:
        self._canvas.bind_all("<MouseWheel>", self._on_wheel)

    def _unbind_wheel(self, _event) -> None:
        self._canvas.unbind_all("<MouseWheel>")

    def _on_wheel(self, event) -> None:
        self._canvas.yview_scroll(int(-event.delta / 120), "units")

    def set_bg(self, bg: str) -> None:
        self._canvas.configure(background=bg)


class Tooltip:
    """Hover tooltip attached to a widget."""

    def __init__(self, widget, text: str, delay_ms: int = 500) -> None:
        self.widget = widget
        self.text = text
        self.delay = delay_ms
        self._after = None
        self._tip = None
        widget.bind("<Enter>", self._schedule)
        widget.bind("<Leave>", self._hide)
        widget.bind("<ButtonPress>", self._hide)

    def _schedule(self, _event=None) -> None:
        self._cancel()
        self._after = self.widget.after(self.delay, self._show)

    def _show(self) -> None:
        if self._tip is not None:
            return
        x = self.widget.winfo_rootx() + 20
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self._tip = tk.Toplevel(self.widget)
        self._tip.wm_overrideredirect(True)
        self._tip.wm_geometry(f"+{x}+{y}")
        tk.Label(self._tip, text=self.text, justify="left",
                 background="#ffffe0", foreground="#000000", relief="solid",
                 borderwidth=1, padx=6, pady=3, wraplength=280).pack()

    def _hide(self, _event=None) -> None:
        self._cancel()
        if self._tip is not None:
            self._tip.destroy()
            self._tip = None

    def _cancel(self) -> None:
        if self._after is not None:
            self.widget.after_cancel(self._after)
            self._after = None
