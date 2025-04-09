# -*- coding: utf-8 -*-
################################################################################
#                                                                              #
#                       Moving Load Beam Analysis                              #                                             #
#                                                                              #
################################################################################

# Version: 1.0 # Initial version
#        - 1.1 # Added UDL support

import sys
import numpy as np
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QLabel, QLineEdit, QPushButton, QTabWidget, QGroupBox, QMessageBox,
                             QRadioButton, QButtonGroup)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QDoubleValidator, QIcon
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import math
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
from typing import Literal


# =============================================================================
# Data Structures / Domain Classes
# =============================================================================

class Beam:
    """Represents the beam being analyzed."""
    def __init__(self, length: float, beam_type: str = "simply_supported"):
        if length <= 0: raise ValueError("Beam length must be positive.")
        self.length = length; self.type = beam_type
    @property
    def L(self): return self.length

class LoadSystemBase(ABC):
    """Abstract Base Class for different load system types."""
    @abstractmethod
    def get_load_type(self) -> str: pass

class PointLoadSystem(LoadSystemBase):
    """Represents the moving two-point load system."""
    LOAD_TYPE = "point_two"
    def __init__(self, w1: float, w2: float, x: float):
        if w1 < 0 or w2 < 0: raise ValueError("Load values cannot be negative.")
        if x < 0: raise ValueError("Distance between loads cannot be negative.")
        self.w1 = w1; self.w2 = w2; self.x = x
    @property
    def W1(self): return self.w1
    @property
    def W2(self): return self.w2
    @property
    def dist_x(self): return self.x
    @property
    def W_total(self): return self.w1 + self.w2
    def get_load_type(self) -> str: return self.LOAD_TYPE

class UDLLoadSystem(LoadSystemBase):
    """Represents a moving Uniformly Distributed Load (UDL)."""
    LOAD_TYPE = "udl"
    def __init__(self, w: float, length_udl: float):
        if w < 0: raise ValueError("UDL intensity (w) cannot be negative.")
        if length_udl <= 0: raise ValueError("UDL length must be positive.")
        self.w = w; self.length = length_udl
    @property
    def intensity(self): return self.w
    @property
    def udl_len(self): return self.length
    def get_load_type(self) -> str: return self.LOAD_TYPE

@dataclass
class AnalysisResults:
    """Holds the results of the beam analysis."""
    load_type: Literal["point_two", "udl", "unknown"] = "unknown"
    Ra_max: float = 0.0; Rb_max: float = 0.0
    SF_max: float = 0.0; y_SF_max: float = 0.0
    BM_01: float | None = 0.0 # Point load specific
    SF_01: float | None = 0.0 # Point load specific
    max_pos_M_mid: float = 0.0 # Max positive Moment @ L/2 (esp. UDL)
    max_pos_SF_mid: float = 0.0 # Max positive Shear @ L/2 (UDL)
    max_neg_SF_mid: float = 0.0 # Max negative Shear @ L/2 (UDL)
    BM_max: float = 0.0 # Absolute max BM (calculated for both types now)
    z_BM_max: float = 0.0 # Location of absolute max BM
    error: str | None = None

# =============================================================================
# Helper Functions for ILD Area Calculation (Unchanged)
# =============================================================================
def _get_ra_ild_ordinate(pos, L): return (L - pos) / L if 0 <= pos <= L and L > 0 else 0.0
def _get_rb_ild_ordinate(pos, L): return pos / L if 0 <= pos <= L and L > 0 else 0.0
def _get_sf_mid_ild_ordinate(pos, L):
    mid = L / 2.0
    if 0 <= pos < mid and L > 0: return -pos / L
    if mid <= pos <= L and L > 0: return (L - pos) / L
    return 0.0
def _get_bm_mid_ild_ordinate(pos, L):
    mid = L / 2.0
    if 0 <= pos <= mid and L > 0: return pos * mid / L
    if mid < pos <= L and L > 0: return (L - pos) * mid / L
    return 0.0
def _calculate_ild_area(ild_ordinate_func, start, end, L):
    if start >= end or L <= 0: return 0.0
    start = max(0, start); end = min(L, end)
    if start >= end : return 0.0
    mid = L / 2.0
    if ild_ordinate_func == _get_sf_mid_ild_ordinate:
        area = 0.0
        if start < mid: eff_end_l=min(end,mid); y_s_l=ild_ordinate_func(start,L); y_e_l=ild_ordinate_func(eff_end_l - 1e-9,L); area += 0.5*(y_s_l+y_e_l)*(eff_end_l-start)
        if end > mid: eff_start_r=max(start,mid); y_s_r=ild_ordinate_func(eff_start_r + 1e-9,L); y_e_r=ild_ordinate_func(end,L); area += 0.5*(y_s_r+y_e_r)*(end-eff_start_r)
        return area
    if ild_ordinate_func == _get_bm_mid_ild_ordinate:
        area = 0.0
        if start <= mid: eff_end_l=min(end,mid); y_s_l=ild_ordinate_func(start,L); y_e_l=ild_ordinate_func(eff_end_l,L); area += 0.5*(y_s_l+y_e_l)*(eff_end_l-start)
        if end > mid: eff_start_r=max(start,mid); y_s_r=ild_ordinate_func(eff_start_r,L); y_e_r=ild_ordinate_func(end,L); area += 0.5*(y_s_r+y_e_r)*(end-eff_start_r)
        return area
    y_start = ild_ordinate_func(start, L); y_end = ild_ordinate_func(end, L)
    area = 0.5 * (y_start + y_end) * (end - start)
    return area

# =============================================================================
# Core Calculation Logic
# =============================================================================

def analyze_beam(beam: Beam, load_system: LoadSystemBase):
    load_type = load_system.get_load_type()
    if beam.type == "simply_supported":
        if load_type == PointLoadSystem.LOAD_TYPE: return _analyze_simply_supported_point(beam, load_system) # type: ignore
        elif load_type == UDLLoadSystem.LOAD_TYPE: return _analyze_simply_supported_udl(beam, load_system) # type: ignore
        else: return AnalysisResults(error=f"Unsupported load type '{load_type}'", load_type='unknown')
    else: return AnalysisResults(error=f"Analysis for beam type '{beam.type}' not implemented.", load_type='unknown')

def _analyze_simply_supported_point(beam: Beam, loads: PointLoadSystem) -> AnalysisResults:
    L = beam.L; W1 = loads.W1; W2 = loads.W2; x = loads.dist_x; W_total = loads.W_total
    res = AnalysisResults(load_type=PointLoadSystem.LOAD_TYPE) # Initializes with 0.0 defaults
    if W_total == 0: return res
    res.Ra_max = W1 + W2 * max(0, L - x) / L
    res.Rb_max = (W1 * max(0, L - x) + W2 * L) / L
    res.SF_max = max(res.Ra_max, res.Rb_max)
    res.y_SF_max = 0.0 if res.Ra_max >= res.Rb_max else L
    if x >= L: res.BM_01 = 0.0
    else: res.BM_01 = (W2 * x / L) * (L - x)
    pos_W1_sf01 = L/2.0; pos_W2_sf01 = pos_W1_sf01 + x
    Ra_for_SF01 = (W1*max(0,L-pos_W1_sf01)+W2*max(0,L-pos_W2_sf01))/L; res.SF_01 = Ra_for_SF01
    # BM_max calculation (unchanged from previous working version)
    if W1==0: res.BM_max=W2*L/4.0; res.z_BM_max=L/2.0
    elif W2==0: res.BM_max=W1*L/4.0; res.z_BM_max=L/2.0
    elif x==0: res.BM_max=W_total*L/4.0; res.z_BM_max=L/2.0
    elif x>=L: res.BM_max=max(W1,W2)*L/4.0; res.z_BM_max=L/2.0
    else:
        R=W_total; a=W2*x/R; z1=L/2.0-a/2.0; M1=-1.0
        if 0<=z1<=L: Ra1=(W1*max(0,L-z1)+W2*max(0,L-(z1+x)))/L; M1=Ra1*z1
        z2=L/2.0+(x-a)/2.0; M2=-1.0
        if 0<=z2<=L: p1_c2=z2-x; Ra2=(W1*max(0,L-p1_c2)+W2*max(0,L-z2))/L; arm_w1=x if p1_c2>=0 else 0; M2=Ra2*z2-W1*arm_w1
        if M1>=0 and M1>=M2: res.BM_max=M1; res.z_BM_max=z1
        elif M2>=0: res.BM_max=M2; res.z_BM_max=z2
        else: res.BM_max=max(W1,W2)*L/4.0; res.z_BM_max=L/2.0
    # Clear UDL specific fields
    res.max_pos_M_mid=None; res.max_pos_SF_mid=None; res.max_neg_SF_mid=None # type: ignore
    return res

def _analyze_simply_supported_udl(beam: Beam, loads: UDLLoadSystem) -> AnalysisResults:
    """ Analysis calculations for UDL, including BM_max. (REVISED) """
    L = beam.L; w = loads.intensity; L_udl = loads.udl_len
    res = AnalysisResults(load_type=UDLLoadSystem.LOAD_TYPE) # Initializes with 0.0 defaults

    # Explicitly set point-load specific fields to None for ALL UDL cases upfront
    res.BM_01 = None
    res.SF_01 = None

    if w == 0: 
        return res

    # --- Max Reactions and Shear ---
    end_ra=min(L_udl,L); area_ra=_calculate_ild_area(_get_ra_ild_ordinate,0,end_ra,L); res.Ra_max=w*area_ra
    start_rb=max(0,L-L_udl); area_rb=_calculate_ild_area(_get_rb_ild_ordinate,start_rb,L,L); res.Rb_max=w*area_rb
    res.SF_max = max(res.Ra_max, res.Rb_max)
    res.y_SF_max = 0.0 if res.Ra_max >= res.Rb_max else L

    # --- Max Effects @ Midspan (using ILDs) ---
    mid=L/2.0
    start_sf_pos=mid; end_sf_pos=min(L,mid+L_udl); area_sf_pos=_calculate_ild_area(_get_sf_mid_ild_ordinate,start_sf_pos,end_sf_pos,L) if start_sf_pos<end_sf_pos else 0.0; res.max_pos_SF_mid = w * area_sf_pos
    end_sf_neg=mid; start_sf_neg=max(0,mid-L_udl); area_sf_neg=_calculate_ild_area(_get_sf_mid_ild_ordinate,start_sf_neg,end_sf_neg,L) if start_sf_neg<end_sf_neg else 0.0; res.max_neg_SF_mid = w * area_sf_neg
    start_m_mid=max(0,mid-L_udl/2.0); end_m_mid=min(L,mid+L_udl/2.0)
    if L_udl >= L: start_m_mid=0; end_m_mid=L
    area_m_mid=_calculate_ild_area(_get_bm_mid_ild_ordinate,start_m_mid,end_m_mid,L) if start_m_mid<end_m_mid else 0.0; res.max_pos_M_mid = w * area_m_mid

    # --- Absolute Max BM for UDL --- ## NEW CALCULATION ##
    if L_udl >= L:
        # Case 1: UDL covers or exceeds the span
        res.BM_max = w * L**2 / 8.0
        res.z_BM_max = L / 2.0
    else:
        # Case 2: UDL is shorter than the span
        # Max moment occurs at L/2 when UDL is centered at L/2
        res.BM_max = (w * L_udl / 8.0) * (2 * L - L_udl)
        res.z_BM_max = L / 2.0

    # --- Clear fields not applicable to UDL analysis ---
    res.BM_01 = None
    res.SF_01 = None

    return res

# =============================================================================
# GUI Application Class (Unchanged)
# =============================================================================
class BeamAnalysisApp(QMainWindow):
    # --- NO CHANGES NEEDED HERE ---
    # --- Keep the entire class definition from the previous version ---
    # --- which includes the UDL input fields and display logic ---
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Osdag® - Beam Analysis (Point Loads & UDL)")
        self.setGeometry(100, 100, 1100, 900)

        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        self.main_layout = QHBoxLayout(self.main_widget)

        self.current_beam: Beam | None = None
        self.current_load_system: LoadSystemBase | None = None
        self.analysis_results: AnalysisResults | None = None

        self.create_input_panel()
        self.create_output_tabs()
        self._update_load_input_visibility() # Set initial visibility

    def create_input_panel(self):
        input_panel = QGroupBox("Input Parameters")
        input_layout = QVBoxLayout()
        validator = QDoubleValidator(); validator.setBottom(0.0)
        pos_validator = QDoubleValidator(0.001, 10000.0, 3)

        beam_input_group = QGroupBox("Beam Properties"); beam_layout = QVBoxLayout()
        self.length_input_field = QLineEdit("10.0"); self.length_input_field.setValidator(pos_validator)
        length_layout = QHBoxLayout(); length_layout.addWidget(QLabel("Beam Length (L) [m]:")); length_layout.addWidget(self.length_input_field); beam_layout.addLayout(length_layout)
        beam_input_group.setLayout(beam_layout); input_layout.addWidget(beam_input_group)

        load_type_group = QGroupBox("Load Type"); load_type_layout = QHBoxLayout()
        self.point_load_radio = QRadioButton("Two Point Loads"); self.udl_radio = QRadioButton("Uniformly Distributed Load (UDL)")
        self.load_type_button_group = QButtonGroup(self); self.load_type_button_group.addButton(self.point_load_radio, 1); self.load_type_button_group.addButton(self.udl_radio, 2)
        self.point_load_radio.setChecked(True); load_type_layout.addWidget(self.point_load_radio); load_type_layout.addWidget(self.udl_radio)
        load_type_group.setLayout(load_type_layout); input_layout.addWidget(load_type_group)
        self.point_load_radio.toggled.connect(self._update_load_input_visibility)

        self.point_load_inputs_group = QGroupBox("Point Load Parameters"); point_load_layout = QVBoxLayout()
        self.w1_input_field = QLineEdit("50.0"); self.w1_input_field.setValidator(validator)
        w1_layout = QHBoxLayout(); w1_layout.addWidget(QLabel("Load W1 [kN]:")); w1_layout.addWidget(self.w1_input_field); point_load_layout.addLayout(w1_layout)
        self.w2_input_field = QLineEdit("30.0"); self.w2_input_field.setValidator(validator)
        w2_layout = QHBoxLayout(); w2_layout.addWidget(QLabel("Load W2 [kN]:")); w2_layout.addWidget(self.w2_input_field); point_load_layout.addLayout(w2_layout)
        self.x_input_field = QLineEdit("2.5"); self.x_input_field.setValidator(validator)
        x_layout = QHBoxLayout(); x_layout.addWidget(QLabel("Distance x [m]:")); x_layout.addWidget(self.x_input_field); point_load_layout.addLayout(x_layout)
        self.point_load_inputs_group.setLayout(point_load_layout); input_layout.addWidget(self.point_load_inputs_group)

        self.udl_inputs_group = QGroupBox("UDL Parameters"); udl_layout = QVBoxLayout()
        self.w_udl_input_field = QLineEdit("10.0"); self.w_udl_input_field.setValidator(validator)
        w_udl_layout = QHBoxLayout(); w_udl_layout.addWidget(QLabel("Intensity w [kN/m]:")); w_udl_layout.addWidget(self.w_udl_input_field); udl_layout.addLayout(w_udl_layout)
        self.len_udl_input_field = QLineEdit("4.0"); self.len_udl_input_field.setValidator(pos_validator)
        len_udl_layout = QHBoxLayout(); len_udl_layout.addWidget(QLabel("UDL Length [m]:")); len_udl_layout.addWidget(self.len_udl_input_field); udl_layout.addLayout(len_udl_layout)
        self.udl_inputs_group.setLayout(udl_layout); input_layout.addWidget(self.udl_inputs_group)

        analyze_btn = QPushButton("Analyze Beam")
        analyze_btn.setStyleSheet("QPushButton { background-color: #4CAF50; color: white; padding: 5px; border-radius: 3px; } QPushButton:hover { background-color: #45a049; }")
        analyze_btn.clicked.connect(self.run_analysis); input_layout.addWidget(analyze_btn)
        input_layout.addStretch(); input_panel.setLayout(input_layout); self.main_layout.addWidget(input_panel, 1)

    def _update_load_input_visibility(self):
        is_point_load = self.point_load_radio.isChecked()
        self.point_load_inputs_group.setVisible(is_point_load)
        self.udl_inputs_group.setVisible(not is_point_load)

    def create_output_tabs(self):
        self.output_tabs = QTabWidget()
        self.results_tab = QWidget(); self.results_layout = QVBoxLayout()
        self.results_display = QLabel("Select load type, enter parameters and click 'Analyze Beam'.")
        self.results_display.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop); self.results_display.setWordWrap(True)
        self.results_layout.addWidget(self.results_display); self.results_tab.setLayout(self.results_layout)
        self.ild_tab = QWidget(); self.ild_layout = QVBoxLayout()
        self.ild_figure = Figure(figsize=(8, 7), dpi=100); self.ild_canvas = FigureCanvas(self.ild_figure)
        self.ild_layout.addWidget(self.ild_canvas); self.ild_tab.setLayout(self.ild_layout)
        self.output_tabs.addTab(self.results_tab, "Summary Results"); self.output_tabs.addTab(self.ild_tab, "Influence Lines (ILDs)")
        self.main_layout.addWidget(self.output_tabs, 3)

    def get_input_and_create_objects(self) -> bool:
        try:
            L = float(self.length_input_field.text()); beam_type = "simply_supported"
            self.current_beam = Beam(length=L, beam_type=beam_type)
            if self.point_load_radio.isChecked():
                W1 = float(self.w1_input_field.text()); W2 = float(self.w2_input_field.text()); x = float(self.x_input_field.text())
                self.current_load_system = PointLoadSystem(w1=W1, w2=W2, x=x)
            elif self.udl_radio.isChecked():
                w = float(self.w_udl_input_field.text()); L_udl = float(self.len_udl_input_field.text())
                self.current_load_system = UDLLoadSystem(w=w, length_udl=L_udl)
            else: raise ValueError("No load type selected.")
            return True
        except ValueError as e: QMessageBox.warning(self, "Input Error", f"Invalid input: {e}"); self.current_beam = None; self.current_load_system = None; return False
        except Exception as e: QMessageBox.critical(self, "Input Error", f"An unexpected error occurred reading inputs: {e}"); self.current_beam = None; self.current_load_system = None; return False

    def run_analysis(self):
        if not self.get_input_and_create_objects(): return
        if self.current_beam is None or self.current_load_system is None: QMessageBox.critical(self, "Internal Error", "Beam or Load System object not created."); return
        self.analysis_results = analyze_beam(self.current_beam, self.current_load_system)
        if self.analysis_results.error:
            QMessageBox.critical(self, "Calculation Error", self.analysis_results.error)
            self.results_display.setText(f"<font color='red'>Error: {self.analysis_results.error}</font>")
            self.ild_figure.clear(); self.ild_canvas.draw(); return
        self.update_results_display()
        self.plot_influence_line_diagrams()

    def update_results_display(self):
        # --- Update this function to display BM_max for UDL ---
        if not self.analysis_results or self.analysis_results.error: return
        res = self.analysis_results; load_type = res.load_type

        results_text = f"<h3>Analysis Results ({load_type.replace('_', ' ').title()})</h3><p>--------------------------------------</p><h4>Maximum Reactions & Shear:</h4>"
        results_text += f"<p><b>Max Reaction at A (Ra_max):</b> {res.Ra_max:.3f} kN</p><p><b>Max Reaction at B (Rb_max):</b> {res.Rb_max:.3f} kN</p><p><b>Absolute Max Shear (SF_max):</b> {res.SF_max:.3f} kN</p><p><b>Location of SF_max (y):</b> {res.y_SF_max:.3f} m from Support A</p><p>--------------------------------------</p>"

        if load_type == PointLoadSystem.LOAD_TYPE:
            ls = self.current_load_system # type: PointLoadSystem
            results_text += "<h4>Point Load Specific Moments & Shear:</h4>"
            if res.BM_01 is not None: results_text += f"<p><b>Moment BM_01:</b> {res.BM_01:.3f} kN·m <br><small><i>(Under W2 when W1@A, x={ls.dist_x:.2f}m)</i></small></p>"
            if res.SF_01 is not None: results_text += f"<p><b>Shear SF_01:</b> {res.SF_01:.3f} kN <br><small><i>(@L/2 when W1@L/2)</i></small></p>"
            # BM_max is now common, displayed later
        elif load_type == UDLLoadSystem.LOAD_TYPE:
            ls = self.current_load_system # type: UDLLoadSystem
            results_text += "<h4>UDL Effects at Midspan (using ILDs):</h4>"
            results_text += f"<p><b>Max Positive Moment at L/2:</b> {res.max_pos_M_mid:.3f} kN·m</p>"
            results_text += f"<p><b>Max Positive Shear at L/2:</b> {res.max_pos_SF_mid:.3f} kN</p>"
            results_text += f"<p><b>Max Negative Shear at L/2:</b> {res.max_neg_SF_mid:.3f} kN</p>"
            results_text += f"<small><i>(Calculated for UDL w={ls.intensity:.2f} kN/m, length={ls.udl_len:.2f}m)</i></small>"

        # Common BM_max display
        results_text += "<p>--------------------------------------</p><h4>Absolute Maximum Bending Moment:</h4>"
        if res.BM_max is not None and res.z_BM_max is not None:
            results_text += f"<p><b>Absolute Max Moment (BM_max):</b> {res.BM_max:.3f} kN·m</p>"
            results_text += f"<p><b>Location of BM_max (z):</b> {res.z_BM_max:.3f} m from Support A</p>"
        else:
             results_text += "<p><i>(Not applicable or calculated for this case)</i></p>"


        results_text += "<p>--------------------------------------</p>"; self.results_display.setText(results_text)

    def plot_influence_line_diagrams(self):
        # --- Plotting logic remains the same (2x2 grid) ---
        self.ild_figure.clear()
        if not self.current_beam or self.current_beam.L <= 0:
            ax = self.ild_figure.add_subplot(111); ax.text(0.5,0.5,'Invalid Beam Data for Plotting', ha='center', va='center', transform=ax.transAxes, fontsize=12, color='red')
            self.ild_canvas.draw(); return
        L = self.current_beam.L; points = np.linspace(0, L, 201)
        ax1 = self.ild_figure.add_subplot(221); ild_ra=(L-points)/L; ax1.plot(points,ild_ra,'b-',label='ILD Ra'); ax1.set_title("ILD Reaction A (Ra)"); ax1.set_ylabel("Value (for P=1)"); ax1.grid(True); ax1.axhline(0, color='k', lw=0.5); ax1.fill_between(points, ild_ra, color='blue', alpha=0.1)
        ax_rb = self.ild_figure.add_subplot(222); ild_rb=points/L; ax_rb.plot(points,ild_rb,'m-',label='ILD Rb'); ax_rb.set_title("ILD Reaction B (Rb)"); ax_rb.set_ylabel("Value (for P=1)"); ax_rb.grid(True); ax_rb.axhline(0, color='k', lw=0.5); ax_rb.fill_between(points, ild_rb, color='magenta', alpha=0.1)
        ax2 = self.ild_figure.add_subplot(223); midspan=L/2.0; ild_sf_corrected=np.zeros_like(points); ild_sf_corrected[points<midspan]=-points[points<midspan]/L; ild_sf_corrected[points>=midspan]=(L-points[points>=midspan])/L # Use >= for right part
        points_sf_plot=np.concatenate(([0],points[points<midspan],[midspan],[midspan],points[points>=midspan],[L])); ild_sf_plot=np.concatenate(([0],-points[points<midspan]/L,[-midspan/L],[(L-midspan)/L],(L-points[points>=midspan])/L,[0] ))
        ax2.plot(points_sf_plot,ild_sf_plot,'r-',label='ILD SF@L/2'); ax2.set_title("ILD Shear @ L/2"); ax2.set_ylabel("Value (for P=1)"); ax2.set_xlabel("Unit Load Position (m)"); ax2.grid(True); ax2.axhline(0, color='k', lw=0.5); ax2.axvline(midspan, color='grey', linestyle='--', linewidth=0.8)
        ax2.fill_between(points_sf_plot,ild_sf_plot, where=ild_sf_plot>=0, color='red', alpha=0.1, interpolate=True); ax2.fill_between(points_sf_plot,ild_sf_plot, where=ild_sf_plot<=0, color='blue', alpha=0.1, interpolate=True)
        ax3 = self.ild_figure.add_subplot(224); ild_bm=np.zeros_like(points); ild_bm[points<=midspan]=(points[points<=midspan]/L)*midspan; ild_bm[points>midspan]=((L-points[points>midspan])/L)*midspan
        ax3.plot(points,ild_bm,'g-',label='ILD BM@L/2'); ax3.set_title("ILD Moment @ L/2"); ax3.set_ylabel("Value (for P=1)"); ax3.set_xlabel("Unit Load Position (m)"); ax3.grid(True); ax3.axhline(0, color='k', lw=0.5); ax3.axvline(midspan, color='grey', linestyle='--', linewidth=0.8)
        ax3.fill_between(points, ild_bm, color='green', alpha=0.1)
        self.ild_figure.tight_layout(pad=2.0); self.ild_canvas.draw()

# =============================================================================
# Main Execution Block
# =============================================================================
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = BeamAnalysisApp()
    window.show()
    sys.exit(app.exec())
# --- END OF FINAL beam.py ---