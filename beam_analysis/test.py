import unittest
import sys
import numpy as np
from PyQt6.QtWidgets import QApplication, QMessageBox
# Import classes and functions
from beam import (Beam, PointLoadSystem, UDLLoadSystem, LoadSystemBase,
                  AnalysisResults, analyze_beam, BeamAnalysisApp,
                   _get_ra_ild_ordinate, _get_rb_ild_ordinate,
                  _get_sf_mid_ild_ordinate, _get_bm_mid_ild_ordinate,
                  _calculate_ild_area )
import io
import math
import os
import datetime
import traceback # To format tracebacks nicely

# --- Configuration ---
RESULTS_FILENAME = "test_results_report.md" # New report filename

# --- Setup QApplication ---
_app_created = False
try:
    app = QApplication.instance();
    if not app: app = QApplication(sys.argv); _app_created = True
except Exception as e: print(f"Warning: QApplication issue ({e}). GUI tests skip."); app = None


# =============================================================================
# Custom Test Result Class to Capture More Info
# =============================================================================
class DetailedTestResult(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.test_details = [] # Store {'name':..., 'status':..., 'details':...}

    def addSuccess(self, test):
        super().addSuccess(test)
        self.test_details.append({'name': test.id(), 'status': 'PASS', 'details': 'OK'})

    def addError(self, test, err):
        super().addError(test, err)
        err_type, err_val, err_tb = err
        tb_str = "".join(traceback.format_exception(err_type, err_val, err_tb))
        self.test_details.append({'name': test.id(), 'status': 'ERROR', 'details': f"{err_type.__name__}: {err_val}\n```\n{tb_str}\n```"})

    def addFailure(self, test, err):
        super().addFailure(test, err)
        err_type, err_val, err_tb = err
        # Get specific assertion message if possible
        # This can be complex, simplifying for now
        details = f"{err_type.__name__}: {err_val}"
        self.test_details.append({'name': test.id(), 'status': 'FAIL', 'details': details})

    def addSkip(self, test, reason):
        super().addSkip(test, reason)
        self.test_details.append({'name': test.id(), 'status': 'SKIP', 'details': reason})

# =============================================================================
# Test Cases for Point Load Calculation Logic
# (Adding parameter logging where feasible)
# =============================================================================
class PointLoadAnalysisTests(unittest.TestCase):

    def _run_test(self, L, W1, W2, x):
        """ Helper to run analysis and store params for reporting """
        self.beam = Beam(length=L)
        self.loads = PointLoadSystem(w1=W1, w2=W2, x=x)
        self.results = analyze_beam(self.beam, self.loads)
        # Store params for potential use in reporting (if test fails)
        self._test_params = f"L={L}, W1={W1}, W2={W2}, x={x}"

    def test_default_case_point(self):
        self._run_test(L=10.0, W1=50.0, W2=30.0, x=2.5)
        self.assertEqual(self.results.load_type, PointLoadSystem.LOAD_TYPE); self.assertIsNone(self.results.error)
        self.assertAlmostEqual(self.results.Ra_max, 72.5, delta=0.01); self.assertAlmostEqual(self.results.Rb_max, 67.5, delta=0.01)
        # ... rest of assertions ...
        self.assertAlmostEqual(self.results.BM_max, 164.2578125, delta=0.01); self.assertAlmostEqual(self.results.z_BM_max, 4.53125, delta=0.01)
        self.assertIsNone(self.results.max_pos_M_mid); self.assertIsNone(self.results.max_pos_SF_mid); self.assertIsNone(self.results.max_neg_SF_mid)

    def test_single_load_w1_point(self):
        self._run_test(L=10.0, W1=50.0, W2=0.0, x=2.5)
        # ... assertions ...
        self.assertAlmostEqual(self.results.Ra_max, 50.0, delta=0.01); self.assertAlmostEqual(self.results.Rb_max, 37.5, delta=0.01)
        self.assertAlmostEqual(self.results.BM_max, 125.0, delta=0.01); self.assertAlmostEqual(self.results.z_BM_max, 5.0, delta=0.01)
        self.assertIsNone(self.results.max_pos_M_mid); self.assertIsNone(self.results.max_pos_SF_mid); self.assertIsNone(self.results.max_neg_SF_mid)

    def test_loads_far_apart_point(self):
        self._run_test(L=10.0, W1=50.0, W2=30.0, x=12.0)
        # ... assertions ...
        self.assertAlmostEqual(self.results.Ra_max, 50.0, delta=0.01); self.assertAlmostEqual(self.results.Rb_max, 30.0, delta=0.01)
        self.assertAlmostEqual(self.results.BM_max, 125.0, delta=0.01); self.assertAlmostEqual(self.results.z_BM_max, 5.0, delta=0.01)
        self.assertIsNone(self.results.max_pos_M_mid); self.assertIsNone(self.results.max_pos_SF_mid); self.assertIsNone(self.results.max_neg_SF_mid)

    def test_point_load_at_midspan(self):
        L=10.0; W1=100; W2=0; x=1.0
        self._run_test(L=L, W1=W1, W2=W2, x=x)
        exp_SF_01 = (W1 * (L/2.0)) / L
        self.assertAlmostEqual(self.results.SF_01, exp_SF_01, delta=0.01)

    def test_coincident_point_loads_robust(self):
        L=10.0; W1=50; W2=30; x=0
        self._run_test(L=L, W1=W1, W2=W2, x=x)
        total_load = W1+W2
        self.assertAlmostEqual(self.results.Ra_max, total_load, delta=0.01)
        self.assertAlmostEqual(self.results.Rb_max, total_load, delta=0.01)
        self.assertAlmostEqual(self.results.SF_01, total_load/2.0, delta=0.01)
        self.assertAlmostEqual(self.results.BM_max, total_load*L/4.0, delta=0.01)

# =============================================================================
# Test Cases for UDL Calculation Logic
# =============================================================================
class UDLAnalysisTests(unittest.TestCase):

    def _run_test(self, L, w, L_udl):
        """ Helper to run analysis and store params for reporting """
        self.beam = Beam(length=L)
        self.loads = UDLLoadSystem(w=w, length_udl=L_udl)
        self.results = analyze_beam(self.beam, self.loads)
        self._test_params = f"L={L}, w={w}, L_udl={L_udl}"

    def test_udl_full_span(self):
        L=10.0; w=5.0; L_udl=12.0
        self._run_test(L=L, w=w, L_udl=L_udl)
        self.assertEqual(self.results.load_type, UDLLoadSystem.LOAD_TYPE); self.assertIsNone(self.results.error)
        exp_reaction=25.0; exp_m_mid=62.5; exp_neg_sf_mid=-6.25; exp_pos_sf_mid=6.25; exp_bm_max=62.5; exp_z_bm_max=5.0
        self.assertAlmostEqual(self.results.Ra_max, exp_reaction, delta=0.01); self.assertAlmostEqual(self.results.Rb_max, exp_reaction, delta=0.01)
        self.assertAlmostEqual(self.results.max_neg_SF_mid, exp_neg_sf_mid, delta=0.01); self.assertAlmostEqual(self.results.BM_max, exp_bm_max, delta=0.01)
        self.assertAlmostEqual(self.results.z_BM_max, exp_z_bm_max, delta=0.01)
        self.assertIsNone(self.results.BM_01); self.assertIsNone(self.results.SF_01)

    def test_udl_short_centered(self):
        L=10.0; w=8.0; L_udl=4.0
        self._run_test(L=L, w=w, L_udl=L_udl)
        self.assertEqual(self.results.load_type, UDLLoadSystem.LOAD_TYPE); self.assertIsNone(self.results.error)
        exp_ra_max=25.6; exp_rb_max=25.6; exp_max_m_mid=64.0; exp_max_pos_sf=9.6; exp_max_neg_sf=-9.6; exp_bm_max=64.0; exp_z_bm_max=5.0
        self.assertAlmostEqual(self.results.Ra_max, exp_ra_max, delta=0.01); self.assertAlmostEqual(self.results.Rb_max, exp_rb_max, delta=0.01)
        self.assertAlmostEqual(self.results.max_neg_SF_mid, exp_max_neg_sf, delta=0.01)
        self.assertAlmostEqual(self.results.BM_max, exp_bm_max, delta=0.01); self.assertAlmostEqual(self.results.z_BM_max, exp_z_bm_max, delta=0.01)
        self.assertIsNone(self.results.BM_01); self.assertIsNone(self.results.SF_01)

    def test_udl_equal_span(self):
        L=8.0; w=10.0; L_udl=8.0
        self._run_test(L=L, w=w, L_udl=L_udl)
        self.assertEqual(self.results.load_type, UDLLoadSystem.LOAD_TYPE); self.assertIsNone(self.results.error)
        exp_reaction=40.0; exp_bm_max=80.0; exp_z_bm_max=4.0
        self.assertAlmostEqual(self.results.Ra_max, exp_reaction, delta=0.01); self.assertAlmostEqual(self.results.Rb_max, exp_reaction, delta=0.01)
        self.assertAlmostEqual(self.results.BM_max, exp_bm_max, delta=0.01); self.assertAlmostEqual(self.results.z_BM_max, exp_z_bm_max, delta=0.01)

    def test_udl_zero_intensity(self):
        L=10.0; w=0.0; L_udl=5.0
        self._run_test(L=L, w=w, L_udl=L_udl)
        self.assertEqual(self.results.load_type, UDLLoadSystem.LOAD_TYPE); self.assertIsNone(self.results.error)
        self.assertAlmostEqual(self.results.Ra_max, 0.0); self.assertAlmostEqual(self.results.Rb_max, 0.0)
        self.assertAlmostEqual(self.results.BM_max, 0.0); self.assertAlmostEqual(self.results.z_BM_max, 0.0)
        self.assertIsNone(self.results.BM_01); self.assertIsNone(self.results.SF_01) # Check these are correctly None

    def test_udl_starts_at_support_a(self):
        L=10.0; w=6.0; L_udl=3.0
        self._run_test(L=L, w=w, L_udl=L_udl)
        expected_ra_max = 15.3
        self.assertAlmostEqual(self.results.Ra_max, expected_ra_max, delta=0.01)

    def test_udl_ends_at_support_b(self):
        L=10.0; w=6.0; L_udl=3.0
        self._run_test(L=L, w=w, L_udl=L_udl)
        expected_rb_max = 15.3
        self.assertAlmostEqual(self.results.Rb_max, expected_rb_max, delta=0.01)

    def test_udl_starts_at_midspan(self):
        L=10.0; w=4.0; L_udl=2.0
        self._run_test(L=L, w=w, L_udl=L_udl)
        expected_max_pos_sf = 3.2
        self.assertAlmostEqual(self.results.max_pos_SF_mid, expected_max_pos_sf, delta=0.01)

    def test_udl_ends_at_midspan(self):
        L=10.0; w=4.0; L_udl=2.0
        self._run_test(L=L, w=w, L_udl=L_udl)
        expected_max_neg_sf = -3.2
        self.assertAlmostEqual(self.results.max_neg_SF_mid, expected_max_neg_sf, delta=0.01)

    def test_udl_very_short(self):
        L=10.0; w=1000.0; L_udl=0.01
        self._run_test(L=L, w=w, L_udl=L_udl)
        self.assertAlmostEqual(self.results.Ra_max, 10.0, delta=0.1)
        self.assertAlmostEqual(self.results.Rb_max, 10.0, delta=0.1)
        self.assertAlmostEqual(self.results.BM_max, 25.0, delta=0.1)

    def test_udl_invalid_inputs(self):
        # Direct checks, no params needed in report for these
        with self.assertRaises(ValueError): Beam(length=0.0)
        with self.assertRaises(ValueError): Beam(length=-5.0)
        with self.assertRaises(ValueError): UDLLoadSystem(w=-2.0, length_udl=5.0)
        with self.assertRaises(ValueError): UDLLoadSystem(w=2.0, length_udl=0.0)
        with self.assertRaises(ValueError): UDLLoadSystem(w=2.0, length_udl=-5.0)

# =============================================================================
# Test Cases for GUI Input Validation (Minimal - No params/results needed)
# =============================================================================
class BeamAnalysisInputTests(unittest.TestCase):
    # ... (Keep class definition as is - no easy way to log GUI params in table) ...
    app_instance = None
    @classmethod
    def setUpClass(cls):
        if app:
            try: print("Creating BeamAnalysisApp instance for GUI tests..."); cls.app_instance = BeamAnalysisApp(); QApplication.processEvents()
            except Exception as e: print(f"\nWarning: Could not create BeamAnalysisApp for GUI tests ({e})."); cls.app_instance = None
        else: print("Skipping GUI test setup: QApplication unavailable.")
    @classmethod
    def tearDownClass(cls):
        if cls.app_instance: print("Closing BeamAnalysisApp instance."); cls.app_instance.close(); QApplication.processEvents(); cls.app_instance = None
    @unittest.skipIf(app is None, "QApplication could not be instantiated")
    def test_gui_valid_inputs_point(self): self.assertIsNotNone(self.app_instance); self.app_instance.point_load_radio.setChecked(True); self.app_instance._update_load_input_visibility(); QApplication.processEvents(); self.app_instance.length_input_field.setText("12.5"); self.app_instance.w1_input_field.setText("100"); self.app_instance.w2_input_field.setText("50"); self.app_instance.x_input_field.setText("3"); self.assertTrue(self.app_instance.get_input_and_create_objects()); self.assertIsInstance(self.app_instance.current_load_system, PointLoadSystem)
    @unittest.skipIf(app is None, "QApplication could not be instantiated")
    def test_gui_valid_inputs_udl(self): self.assertIsNotNone(self.app_instance); self.app_instance.udl_radio.setChecked(True); self.app_instance._update_load_input_visibility(); QApplication.processEvents(); self.app_instance.length_input_field.setText("15.0"); self.app_instance.w_udl_input_field.setText("12"); self.app_instance.len_udl_input_field.setText("6.5"); self.assertTrue(self.app_instance.get_input_and_create_objects()); self.assertIsInstance(self.app_instance.current_load_system, UDLLoadSystem)
    @unittest.skipIf(app is None, "QApplication could not be instantiated")
    def test_gui_invalid_length(self): self.assertIsNotNone(self.app_instance); self.app_instance.length_input_field.setText("0"); self.assertFalse(self.app_instance.get_input_and_create_objects())
    @unittest.skipIf(app is None, "QApplication could not be instantiated")
    def test_gui_non_numeric_input_point(self): self.assertIsNotNone(self.app_instance); self.app_instance.point_load_radio.setChecked(True); self.app_instance._update_load_input_visibility(); QApplication.processEvents(); self.app_instance.w1_input_field.setText("abc"); self.assertFalse(self.app_instance.get_input_and_create_objects())
    @unittest.skipIf(app is None, "QApplication could not be instantiated")
    def test_gui_invalid_udl_length(self): self.assertIsNotNone(self.app_instance); self.app_instance.udl_radio.setChecked(True); self.app_instance._update_load_input_visibility(); QApplication.processEvents(); self.app_instance.len_udl_input_field.setText("0"); self.assertFalse(self.app_instance.get_input_and_create_objects())


# =============================================================================
# Main Test Execution and Reporting
# =============================================================================
if __name__ == "__main__":
    print("Running robust tests...")

    # Create test suites
    suite_point = unittest.TestLoader().loadTestsFromTestCase(PointLoadAnalysisTests)
    suite_udl = unittest.TestLoader().loadTestsFromTestCase(UDLAnalysisTests)
    suite_gui = unittest.TestLoader().loadTestsFromTestCase(BeamAnalysisInputTests)
    full_suite = unittest.TestSuite([suite_point, suite_udl, suite_gui])

    # Use the custom result class
    runner = unittest.TextTestRunner(resultclass=DetailedTestResult, stream=sys.stderr, verbosity=2) # Print runner output to stderr

    print("-" * 70)
    result = runner.run(full_suite) # Run tests and collect detailed results
    print("-" * 70)


    # --- Generate Markdown Report ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    report_path = os.path.join(script_dir, RESULTS_FILENAME)
    print(f"Generating test report: {report_path}")

    total_run = result.testsRun
    failures = len(result.failures)
    errors = len(result.errors)
    skipped = len(result.skipped)
    passed = total_run - failures - errors - skipped
    success = failures == 0 and errors == 0

    try:
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(f"# Beam Analysis Test Report\n\n")
            f.write(f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            f.write(f"## Summary\n\n")
            f.write("| Status      | Count |\n")
            f.write("| :---------- | ----: |\n")
            f.write(f"| **PASSED**  | {passed} |\n")
            f.write(f"| **FAILED**  | {failures} |\n")
            f.write(f"| **ERRORS**  | {errors} |\n")
            f.write(f"| **SKIPPED** | {skipped} |\n")
            f.write(f"| **TOTAL**   | {total_run} |\n")
            f.write(f"\n**Overall Status:** {'✅ PASS' if success else '❌ FAIL'}\n\n")


            f.write(f"## Detailed Test Results\n\n")
            f.write("| Test Case                                                    | Status | Details / Parameters |\n")
            f.write("| :----------------------------------------------------------- | :----- | :------------------- |\n")

            # Access captured details from the custom result object
            for detail in result.test_details:
                name = detail['name']
                status = detail['status']
                details_text = detail['details']
                # Try to get params if available (simple approach)
                params_text = ""
                # This relies on the test method having stored params in self._test_params
                # It's a bit hacky - requires finding the test object corresponding to the name
                # For simplicity, only show params on failure/error for calculation tests
                if status in ['FAIL', 'ERROR'] and ('PointLoadAnalysisTests' in name or 'UDLAnalysisTests' in name):
                    # Find the test method instance (if possible - complex in standard runner)
                    # This part is difficult without deeper integration or storing params directly
                    # in the result object, which TextTestResult doesn't easily support.
                    # We'll just show the error detail here. A more advanced framework
                    # like pytest might make parameter reporting easier.
                     pass # Can't easily get params back here

                # Escape pipe characters in details for Markdown table
                details_text = details_text.replace('|', '\\|').replace('\n', '<br>')
                # Truncate long details?
                max_detail_len = 200
                if len(details_text) > max_detail_len:
                    details_text = details_text[:max_detail_len] + "..."

                status_icon = {'PASS': '✅', 'FAIL': '❌', 'ERROR': '🔥', 'SKIP': '⏭️'}.get(status, '❓')
                f.write(f"| `{name}` | {status_icon} {status} | {details_text} |\n")

        print("Report generated successfully.")

    except IOError as e:
        print(f"Error writing report file: {e}")

    exit_code = 0 if success else 1
    if _app_created and app: print("Quitting created QApplication instance."); app.exit()
    sys.exit(exit_code)