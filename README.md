# Moving Load Beam Analysis using Influence Lines

## Overview

This project provides a Python application with a graphical user interface (GUI) to analyze a simply supported beam subjected to a moving two-point load system. It calculates key structural responses like maximum reactions, shear forces, and bending moments using principles derived from Influence Line Diagrams (ILDs). The application is built using PyQt6 for the GUI and Matplotlib for plotting ILDs.

The primary goal is to demonstrate the calculation of critical design values for beams under moving loads and to illustrate the concept of Influence Lines.

**Target Audience:** Civil/Structural Engineering students, engineers, or anyone interested in structural analysis and programming.

## Features

*   Calculates structural responses for a simply supported beam under two moving point loads (W1, W2) separated by a fixed distance (x).
*   Determines:
    *   Maximum Reaction at Support A (Ra_max)
    *   Maximum Reaction at Support B (Rb_max)
    *   Shear Force at mid-span (L/2) when W1 is positioned at mid-span (SF\_01)
    *   Bending Moment under W2 when W1 is positioned at Support A (BM\_01)
    *   Absolute Maximum Shear Force (SF\_max) and its location (y)
    *   Absolute Maximum Bending Moment (BM\_max) and its location (z)
*   Plots Influence Line Diagrams (ILDs) for:
    *   Reaction at Support A (Ra)
    *   Shear Force at Mid-span (V @ L/2)
    *   Bending Moment at Mid-span (M @ L/2)
*   User-friendly GUI for inputting parameters (Beam Length L, Loads W1 & W2, Distance x).
*   Clear display of results and ILD plots in separate tabs.
*   Includes unit tests to verify calculation logic.

## Screenshots

*(Add screenshots of your application here)*

**Input Panel & Results Tab:**
`[Screenshot showing the input fields and the results tab with calculated values]`

**Influence Lines Tab:**
`[Screenshot showing the plotted ILDs for Ra, SF@L/2, and BM@L/2]`

## Core Concepts: Influence Line Diagrams (ILDs) for Moving Loads

Understanding how structures respond to loads that change position (like trucks on a bridge, cranes on a beam) is crucial in structural design. We need to find the *maximum* possible force or moment that a specific point in the structure will experience as the load traverses it. Calculating this for every possible load position is tedious. This is where Influence Line Diagrams become powerful.

### What is an Influence Line Diagram (ILD)?

An Influence Line Diagram (ILD) is a graph that shows the variation of a *specific function* (like reaction at support A, shear force at point C, or bending moment at point D) *at a fixed point* on the structure as a **single unit load (value = 1)** moves across the entire length of the structure.

*   **X-axis:** Position of the moving unit load.
*   **Y-axis:** Value of the function (Reaction, Shear, Moment) *at the specific point of interest* caused by the unit load at position X.

**Key Idea:** The ILD tells you how influential the unit load's position is on the specific response you care about (e.g., how much does placing a 1 kN load at `x=L/4` contribute to the reaction at support A?).

### Why Use ILDs for Moving Loads?

Imagine multiple loads (like the wheels of a truck, represented here by W1 and W2) moving together. To find the maximum reaction at A, you'd theoretically have to:
1.  Place the load system (W1, W2) at position 1. Calculate Ra.
2.  Move the system slightly to position 2. Calculate Ra.
3.  Repeat for *all* possible positions.
4.  Find the maximum Ra from all calculations.

This is inefficient. ILDs provide a visual and mathematical shortcut.

### How ILDs Work: Using the Diagram

The ordinate (y-value) of an ILD at a specific position `x` represents the value of the function (e.g., Ra) when a unit load is placed at that position `x`.

For a system of multiple loads (P1, P2, P3...) acting simultaneously, the **total value of the function** is found by summing the contribution of each load, using the principle of superposition:

`Total Effect = P1 * y1 + P2 * y2 + P3 * y3 + ... = Σ (Pi * yi)`

Where:
*   `Pi` is the magnitude of the i-th load.
*   `yi` is the ordinate (y-value) of the ILD directly *under* the position of load `Pi`.

### Finding Maximum Effects using ILDs

To find the maximum value of a function (e.g., max Ra, max positive moment at C) caused by a moving load system (like W1, W2):

1.  **Identify the ILD:** Draw or obtain the ILD for the specific function you are interested in (e.g., ILD for Ra).
2.  **Position the Loads:** Place the moving load system (W1, W2) onto the ILD diagram. Mentally (or algorithmically) slide the load system across the ILD.
3.  **Maximize the Sum:** Find the position of the load system that maximizes the sum `Σ (Pi * yi)`.
    *   For maximum *positive* effects (like max Ra, max positive BM): Position the loads over the largest positive ordinates of the ILD. Often, this involves placing the heaviest load(s) near the peak(s) of the ILD.
    *   For maximum *negative* effects (like max negative BM, max negative shear): Position loads over the largest negative ordinates.
    *   For maximum *absolute* shear: Check both maximum positive and maximum negative possibilities.

### Specific ILDs and Calculations in this Project

*   **ILD for Reaction at A (Ra):**
    *   *Shape:* Triangle, peak value of 1 at Support A (x=0), value of 0 at Support B (x=L). Equation: `y = (L-x)/L`.
    *   *Max Ra:* To maximize `Ra = W1*y1 + W2*y2`, we need the largest ordinates under W1 and W2. This clearly occurs when W1 is placed directly over Support A (y1=1), and W2 is at distance x (y2 = (L-x)/L). The code calculates `Ra_max` using this principle directly: `Ra_max = W1 * 1 + W2 * (L - x) / L`.
*   **ILD for Shear at Midspan (V @ L/2):**
    *   *Shape:* Two triangles. From x=0 to x=L/2, it goes linearly from 0 to -0.5 (or -L/(2L) = -0.5). From x=L/2 to x=L, it goes linearly from +0.5 to 0. There's a sharp jump at x=L/2.
    *   *SF_01 Calculation:* The requirement `SF_01` asks for the shear *at* L/2 when W1 is placed *exactly* at L/2. The code calculates the reaction Ra for this specific load position (`Ra = (W1*(L-L/2) + W2*max(0, L-(L/2+x))) / L`) and reports this Ra as the shear just to the left of the loads at midspan. This aligns with using the ILD principle for a specific load configuration.
*   **ILD for Moment at Midspan (M @ L/2):**
    *   *Shape:* Triangle/Parabola, peak value of L/4 at the mid-span (x=L/2), value of 0 at supports. Equation: `y = x/2` for `x <= L/2` and `y = (L-x)/2` for `x > L/2`.
    *   *Max M @ L/2:* To maximize moment *at this specific point*, loads W1 and W2 would be placed near the peak (L/2) to capture the largest ordinates. The code *plots* this ILD but calculates the *absolute* maximum moment (BM_max) separately, which might occur at a different location and under a different load.
*   **Absolute Maximum Shear (SF_max):** For a simply supported beam, the absolute maximum shear occurs when one of the heaviest loads is placed infinitesimally close to a support, maximizing the reaction at that support. The ILD for reactions confirms this (peak value of 1 at the support). Therefore, `SF_max` is simply the larger of `Ra_max` and `Rb_max`. The location `y` is 0 (Support A) if `Ra_max` governs, or L (Support B) if `Rb_max` governs.
*   **Absolute Maximum Bending Moment (BM_max):** Finding the absolute maximum bending moment under a *series* of moving loads is more complex. It generally occurs *under* one of the loads (usually the one closer to the resultant of the load group) when the load system is positioned such that the **centerline of the beam bisects the distance between that specific load and the resultant of the load group**. The code implements formulas derived from this principle in the `analyze_beam_calculations` function to find `BM_max` and its location `z`. This calculation uses the underlying principles demonstrated by ILDs but doesn't require plotting the ILD for that specific (and initially unknown) critical moment location. The `BM_01` calculation (moment under W2 when W1 is at A) is a specific case analysis, not the absolute maximum.

## Installation

1.  **Clone the repository:**
    ```bash
    git clone <your-repository-url>
    cd <repository-directory>
    ```
2.  **Create and activate a virtual environment (recommended):**
    ```bash
    # For Linux/macOS
    python3 -m venv venv
    source venv/bin/activate

    # For Windows
    python -m venv venv
    .\venv\Scripts\activate
    ```
3.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

## Usage

1.  Make sure your virtual environment is activated.
2.  Run the main application script:
    ```bash
    python beam.py
    ```
3.  The GUI window will appear.
4.  Enter the required parameters in the "Input Parameters" section:
    *   Beam Length (L) in meters.
    *   Moving Load W1 in kN.
    *   Moving Load W2 in kN.
    *   Distance W1-W2 (x) in meters.
5.  Click the "Analyze Beam" button.
6.  View the calculated results in the "Summary Results" tab.
7.  View the plotted Influence Line Diagrams in the "Influence Lines (ILDs)" tab.

## File Structure
```
└── beam_analysis
    ├── beam.py
    ├── test.py 
    ├── test_results.md
    └── README.md 

```