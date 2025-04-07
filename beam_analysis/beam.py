

import sys
import numpy as np
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QLabel, QLineEdit, QPushButton, QTabWidget, QGroupBox, QMessageBox)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QDoubleValidator, QIcon 
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import math

# =============================================================================
#                            Core Calculation Logic
# =============================================================================

def analyze_beam_calculations(L, W1, W2, x):
    """
    Calculates reactions, shear forces, and bending moments for a simply supported beam
    subjected to a moving two-point load (W1, W2 separated by distance x).

    Args:
        L (float): Length of the beam (m). Must be > 0.
        W1 (float): Magnitude of the first moving load (kN). Must be >= 0.
        W2 (float): Magnitude of the second moving load (kN). Must be >= 0.
        x (float): Distance between W1 and W2 (m). Must be >= 0.

    Returns:
        dict: A dictionary containing the calculated results:
              {'Ra_max', 'Rb_max', 'BM_01', 'SF_01', 'SF_max', 'y_SF_max',
               'BM_max', 'z_BM_max', 'error': None}
              or {'error': 'Error message string'} if inputs are invalid.
    """
    results = {}

    # --- Input Validation ---
    if L <= 0:
        # Use the exact error string expected by the test
        return {'error': "Beam length (L) must be positive."}
    if W1 < 0 or W2 < 0:
        return {'error': "Load values (W1, W2) cannot be negative."}
    if x < 0:
        return {'error': "Distance between loads (x) cannot be negative."}
    if W1 == 0 and W2 == 0:
        # Handle zero load case gracefully
        return {
            'Ra_max': 0.0, 'Rb_max': 0.0, 'BM_01': 0.0, 'SF_01': 0.0,
            'SF_max': 0.0, 'y_SF_max': 0.0, 'BM_max': 0.0, 'z_BM_max': 0.0,
            'error': None
        }

    # --- Maximum Reactions ---
    # Max Ra: Occurs when W1 is at support A (position 0). W2 is at x.
    # Ra = W1*L/L + W2*max(0, L-x)/L
    results['Ra_max'] = W1 + W2 * max(0, L - x) / L # L>0 check already done

    # Max Rb: Occurs when W2 is at support B (position L). W1 is at L-x.
    # Rb = W1*max(0, L-x)/L + W2*L/L  <- CORRECTED FORMULA
    results['Rb_max'] = (W1 * max(0, L - x) + W2 * L) / L


    # --- BM_01: Bending moment when W1 is at support A (0 m) ---
    # Moment is calculated under W2, which is at position x.
    if x >= L: # W2 is off the beam
         results['BM_01'] = 0.0
    else:
        # Reaction Rb when W1 is at 0: Rb = (W1 * 0 + W2 * x) / L = W2 * x / L
        # Moment under W2 = Rb * (L - x)
        results['BM_01'] = (W2 * x / L) * (L - x)

    # --- Shear Force Calculations ---
    # SF_max: Absolute maximum shear occurs at a support and equals the max reaction
    results['SF_max'] = max(results['Ra_max'], results['Rb_max'])
    results['y_SF_max'] = 0.0 if results['Ra_max'] >= results['Rb_max'] else L # Location of SF_max

    # SF_01: Shear force at midspan (L/2)
    # Calculated for the specific position where W1 is AT L/2
    # L > 0 check already performed at the start
    pos_W1_sf01 = L / 2.0
    pos_W2_sf01 = pos_W1_sf01 + x
    # Ra for this specific position:
    Ra_for_SF01 = (W1 * max(0, L - pos_W1_sf01) + W2 * max(0, L - pos_W2_sf01)) / L
    # SF at L/2 (just to the left of W1 if it's exactly at L/2) is Ra
    results['SF_01'] = Ra_for_SF01


    # --- Absolute Maximum Bending Moment (BM_max) ---
    W_total = W1 + W2
    # L > 0 and W_total > 0 checks already done or handled
    if W1 == 0: # Only W2 acts
        results['BM_max'] = W2 * L / 4.0
        results['z_BM_max'] = L / 2.0
    elif W2 == 0: # Only W1 acts
        results['BM_max'] = W1 * L / 4.0
        results['z_BM_max'] = L / 2.0
    elif x == 0: # Loads are coincident
        results['BM_max'] = (W1 + W2) * L / 4.0
        results['z_BM_max'] = L / 2.0
    elif x >= L: # Loads cannot act together for max moment near center
         results['BM_max'] = max(W1, W2) * L / 4.0 # Governed by heavier single load
         results['z_BM_max'] = L / 2.0
    else:
        # Two distinct loads, 0 < x < L
        # Find resultant R and its distance 'a' from W1
        R = W1 + W2
        a = W2 * x / R # Distance from W1 towards W2

        # Case 1: Check moment under W1
        # Position W1 such that centerline (L/2) is midway between W1 and R.
        z1 = L / 2.0 - a / 2.0
        M1 = -1.0 # Initialize to invalid moment

        # Check if W1 is on the beam for this position
        if 0 <= z1 <= L:
             pos_W1_case1 = z1
             pos_W2_case1 = z1 + x
             # Calculate Ra for this position
             Ra1 = (W1 * max(0, L - pos_W1_case1) + W2 * max(0, L - pos_W2_case1)) / L
             # Moment under W1
             M1 = Ra1 * z1


        # Case 2: Check moment under W2
        # Position W2 such that centerline (L/2) is midway between W2 and R.
        z2 = L / 2.0 + (x - a) / 2.0
        M2 = -1.0 # Initialize to invalid moment

         # Check if W2 is on the beam for this position
        if 0 <= z2 <= L:
            pos_W2_case2 = z2
            pos_W1_case2 = z2 - x
            # Calculate Ra for this position
            Ra2 = (W1 * max(0, L - pos_W1_case2) + W2 * max(0, L - pos_W2_case2)) / L
            # Moment under W2 = Ra * z2 - W1 * (distance between W1 and W2, if W1 is to the left)
            moment_arm_w1 = x if pos_W1_case2 >= 0 else 0 # Check if W1 is on beam
            M2 = Ra2 * z2 - W1 * moment_arm_w1


        # Determine the absolute maximum moment and its location
        # Prioritize valid moments. If both valid, take the larger one.
        if M1 >= 0 and M1 >= M2: # Handles M2 being invalid (-1) or smaller
            results['BM_max'] = M1
            results['z_BM_max'] = z1
        elif M2 >= 0: # M2 is valid and greater than M1 (or M1 was invalid)
             results['BM_max'] = M2
             results['z_BM_max'] = z2
        else:
             # Fallback if somehow both M1 and M2 are invalid (e.g. < 0)
             # This case is highly unlikely with valid positive loads and L>0, x<L
             # Default to heavier single load case as a safety measure
             results['BM_max'] = max(W1, W2) * L / 4.0
             results['z_BM_max'] = L / 2.0


    # Check for NaN results which might occur in edge cases before returning
    final_results = {}
    final_results['error'] = None # Assume success initially
    for key, value in results.items():
        # Check if value is a float and is NaN
        if isinstance(value, float) and math.isnan(value):
             # Return an error if NaN found
             return {'error': f"Calculation resulted in NaN for {key}."}
        final_results[key] = value # Copy valid results

    return final_results

# =============================================================================
#                             GUI Application Class 
# =============================================================================

class BeamAnalysisApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Osdag® - Beam Analysis with ILD (Moving Load)")
        self.setGeometry(100, 100, 1100, 700) # Adjusted size slightly
        # Consider adding an application icon
        # self.setWindowIcon(QIcon("path/to/your/icon.png"))

        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        self.main_layout = QHBoxLayout(self.main_widget)

        # Initialize variables
        self.beam_length = 0.0
        self.w1 = 0.0
        self.w2 = 0.0
        self.x_distance = 0.0
        self.results = {} # Store results from analyze_beam_calculations

        self.create_input_panel()
        self.create_output_tabs()

    def create_input_panel(self):
        input_panel = QGroupBox("Input Parameters")
        input_layout = QVBoxLayout()

        # Use QDoubleValidator for better input control
        validator = QDoubleValidator()
        validator.setBottom(0.0) # Length, loads, distance >= 0

        # Beam Length
        self.length_input_field = QLineEdit("10.0")
        self.length_input_field.setValidator(validator)
        length_layout = QHBoxLayout()
        length_layout.addWidget(QLabel("Beam Length (L) [m]:"))
        length_layout.addWidget(self.length_input_field)
        input_layout.addLayout(length_layout)

        # Load W1
        self.w1_input_field = QLineEdit("50.0")
        self.w1_input_field.setValidator(validator)
        w1_layout = QHBoxLayout()
        w1_layout.addWidget(QLabel("Moving Load W1 [kN]:"))
        w1_layout.addWidget(self.w1_input_field)
        input_layout.addLayout(w1_layout)

        # Load W2
        self.w2_input_field = QLineEdit("30.0")
        self.w2_input_field.setValidator(validator)
        w2_layout = QHBoxLayout()
        w2_layout.addWidget(QLabel("Moving Load W2 [kN]:"))
        w2_layout.addWidget(self.w2_input_field)
        input_layout.addLayout(w2_layout)

        # Distance between loads
        self.x_input_field = QLineEdit("2.5")
        self.x_input_field.setValidator(validator) # Allow x=0
        x_layout = QHBoxLayout()
        x_layout.addWidget(QLabel("Distance W1-W2 (x) [m]:"))
        x_layout.addWidget(self.x_input_field)
        input_layout.addLayout(x_layout)

        # Analyze button
        analyze_btn = QPushButton("Analyze Beam")
        analyze_btn.setStyleSheet("QPushButton { background-color: #4CAF50; color: white; padding: 5px; border-radius: 3px; } QPushButton:hover { background-color: #45a049; }")
        analyze_btn.clicked.connect(self.run_analysis)
        input_layout.addWidget(analyze_btn)

        input_layout.addStretch()
        input_panel.setLayout(input_layout)
        self.main_layout.addWidget(input_panel, 1) # Input panel takes less space

    def create_output_tabs(self):
        self.output_tabs = QTabWidget()

        # --- Results Tab ---
        self.results_tab = QWidget()
        self.results_layout = QVBoxLayout()
        self.results_display = QLabel("Enter parameters and click 'Analyze Beam'.")
        self.results_display.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.results_display.setWordWrap(True) # Allow text wrapping
        self.results_layout.addWidget(self.results_display)
        self.results_tab.setLayout(self.results_layout)

        # --- ILD Tab ---
        self.ild_tab = QWidget()
        self.ild_layout = QVBoxLayout()
        self.ild_figure = Figure(figsize=(8, 6), dpi=100) # Adjusted figure size
        self.ild_canvas = FigureCanvas(self.ild_figure)
        self.ild_layout.addWidget(self.ild_canvas)
        self.ild_tab.setLayout(self.ild_layout)

        # Add tabs to the widget
        self.output_tabs.addTab(self.results_tab, "Summary Results")
        self.output_tabs.addTab(self.ild_tab, "Influence Lines (ILDs)")
        # Removed the misleading SF/BM Diagram Tab

        self.main_layout.addWidget(self.output_tabs, 3) # Output tabs take more space

    def get_input_values(self):
        """Retrieves and validates input values from GUI fields."""
        try:
            L = float(self.length_input_field.text())
            W1 = float(self.w1_input_field.text())
            W2 = float(self.w2_input_field.text())
            x = float(self.x_input_field.text())

            # Basic validation (more detailed validation in analyze_beam_calculations)
            if L <= 0:
                raise ValueError("Beam length must be positive.")
            # Allow x >= L, handled in calculation logic
            # Allow W1=0 or W2=0

            self.beam_length = L
            self.w1 = W1
            self.w2 = W2
            self.x_distance = x
            return True

        except ValueError as e:
            QMessageBox.warning(self, "Input Error", f"Invalid input: {e}\nPlease enter valid numbers.")
            return False
        except Exception as e: # Catch other potential errors
             QMessageBox.critical(self, "Input Error", f"An unexpected error occurred reading inputs: {e}")
             return False

    def run_analysis(self):
        """Gets inputs, runs calculations, and updates outputs."""
        if not self.get_input_values():
            return # Stop if inputs are invalid

        # Perform calculations using the refactored function
        self.results = analyze_beam_calculations(
            self.beam_length, self.w1, self.w2, self.x_distance
        )

        # Check if calculation returned an error
        if self.results.get('error'):
            QMessageBox.critical(self, "Calculation Error", self.results['error'])
            # Clear previous results display maybe?
            self.results_display.setText(f"<font color='red'>Error: {self.results['error']}</font>")
            # Clear plots or show placeholder?
            self.ild_figure.clear()
            self.ild_canvas.draw()
            return

        # Update the GUI with results
        self.update_results_display()
        self.plot_influence_line_diagrams()
        # Removed call to plot_sf_bm_diagrams

    def update_results_display(self):
        """Formats and displays the calculated results in the Results Tab."""
        if not self.results or self.results.get('error'):
            return # Don't display if no results or error occurred

        # Format results nicely using HTML for basic styling
        results_text = f"""
        <h3>Analysis Results Summary</h3>
        <p>--------------------------------------</p>
        <h4>Maximum Reactions:</h4>
        <p><b>Max Reaction at A (Ra_max):</b> {self.results['Ra_max']:.3f} kN</p>
        <p><b>Max Reaction at B (Rb_max):</b> {self.results['Rb_max']:.3f} kN</p>
        <p>--------------------------------------</p>
        <h4>Shear Forces:</h4>
        <p><b>Shear Force at L/2 (SF_01):</b> {self.results['SF_01']:.3f} kN <br>
           <small><i>(Calculated with W1 positioned exactly at L/2)</i></small></p>
        <p><b>Absolute Max Shear (SF_max):</b> {self.results['SF_max']:.3f} kN</p>
        <p><b>Location of SF_max (y):</b> {self.results['y_SF_max']:.3f} m from Support A</p>
        <p>--------------------------------------</p>
        <h4>Bending Moments:</h4>
        <p><b>Moment BM_01:</b> {self.results['BM_01']:.3f} kN·m <br>
           <small><i>(Calculated under W2 when W1 is at Support A, x={self.x_distance}m)</i></small></p>
        <p><b>Absolute Max Moment (BM_max):</b> {self.results['BM_max']:.3f} kN·m</p>
        <p><b>Location of BM_max (z):</b> {self.results['z_BM_max']:.3f} m from Support A</p>
        <p>--------------------------------------</p>
        """
        self.results_display.setText(results_text)

    def plot_influence_line_diagrams(self):
        """Plots standard Influence Line Diagrams for RA, SF@midspan, BM@midspan."""
        self.ild_figure.clear()
        if self.beam_length <= 0:
            # Add a placeholder message if length is invalid
            ax = self.ild_figure.add_subplot(111)
            ax.text(0.5, 0.5, 'Invalid Beam Length for Plotting',
                    horizontalalignment='center', verticalalignment='center',
                    transform=ax.transAxes, fontsize=12, color='red')
            self.ild_canvas.draw()
            return

        L = self.beam_length
        points = np.linspace(0, L, 201) # More points for smoother curves

        # --- ILD for Reaction at A ---
        ax1 = self.ild_figure.add_subplot(311)
        ild_ra = (L - points) / L if L > 0 else np.zeros_like(points) # Handle L=0 case
        ax1.plot(points, ild_ra, 'b-', label='ILD for Ra')
        ax1.set_title("Influence Line for Reaction at A (Ra)")
        ax1.set_ylabel("Ra value (for 1 unit load)")
        ax1.grid(True)
        ax1.axhline(0, color='black', linewidth=0.5)
        ax1.fill_between(points, ild_ra, color='blue', alpha=0.1) # Fill area

        # --- ILD for Shear Force at Midspan (L/2) ---
        ax2 = self.ild_figure.add_subplot(312)
        midspan = L / 2.0
        # Define points slightly left and right of midspan for clear vertical line
        points_sf = np.linspace(0, L, 201)
        # Standard ILD for SF at section 'c': value is Ra when load is right, -Rb when load is left.
        # At midspan (c=L/2): Ra = (L-x)/L, Rb = x/L
        ild_sf_corrected = np.zeros_like(points_sf)
        if L > 0: # Avoid division by zero
            ild_sf_corrected[points_sf < midspan] = -points_sf[points_sf < midspan] / L # Value is -Rb
            ild_sf_corrected[points_sf > midspan] = (L - points_sf[points_sf > midspan]) / L # Value is Ra
            # Add discontinuity point explicitly for plotting
            points_sf_plot = np.concatenate(([0], points_sf[points_sf < midspan], [midspan], [midspan], points_sf[points_sf > midspan], [L]))
            ild_sf_plot = np.concatenate(([0], -points_sf[points_sf < midspan]/L, [-midspan/L], [(L-midspan)/L], (L - points_sf[points_sf > midspan])/L, [0] ))
        else: # Handle L=0 case
            points_sf_plot = np.array([0])
            ild_sf_plot = np.array([0])


        ax2.plot(points_sf_plot, ild_sf_plot, 'r-', label='ILD for SF @ L/2')
        ax2.set_title("Influence Line for Shear Force at Midspan (SF @ L/2)")
        ax2.set_ylabel("SF value (for 1 unit load)")
        ax2.grid(True)
        ax2.axhline(0, color='black', linewidth=0.5)
        ax2.axvline(midspan, color='grey', linestyle='--', linewidth=0.8)
        ax2.fill_between(points_sf_plot, ild_sf_plot, where=ild_sf_plot>=0, color='red', alpha=0.1, interpolate=True)
        ax2.fill_between(points_sf_plot, ild_sf_plot, where=ild_sf_plot<=0, color='blue', alpha=0.1, interpolate=True)


        # --- ILD for Bending Moment at Midspan (L/2) ---
        ax3 = self.ild_figure.add_subplot(313)
        midspan = L / 2.0
        # BM at c due to unit load at x: Ra*c for x>c, Rb*(L-c) for x<c
        # At c=L/2: Ra = (L-x)/L, Rb = x/L
        ild_bm = np.zeros_like(points)
        if L > 0: # Avoid division by zero
            ild_bm[points <= midspan] = (points[points <= midspan] / L) * (L - midspan) # Rb * (L-c) = x/L * L/2 = x/2
            ild_bm[points > midspan] = ((L - points[points > midspan]) / L) * midspan # Ra * c = (L-x)/L * L/2 = (L-x)/2
            # Peak value is (L/2)/2 = L/4

        ax3.plot(points, ild_bm, 'g-', label='ILD for BM @ L/2')
        ax3.set_title("Influence Line for Bending Moment at Midspan (BM @ L/2)")
        ax3.set_xlabel("Position of Unit Load along Beam (m)")
        ax3.set_ylabel("BM value (for 1 unit load)")
        ax3.grid(True)
        ax3.axhline(0, color='black', linewidth=0.5)
        ax3.axvline(midspan, color='grey', linestyle='--', linewidth=0.8)
        ax3.fill_between(points, ild_bm, color='green', alpha=0.1) # Fill area

        # Adjust layout and draw
        self.ild_figure.tight_layout(pad=2.0) # Add padding
        self.ild_canvas.draw()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = BeamAnalysisApp()
    window.show()
    sys.exit(app.exec())
