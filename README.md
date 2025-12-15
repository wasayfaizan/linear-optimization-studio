<div align="center">

<h1>🔷 Linear Optimization Studio</h1>
<p><i>A modern, interactive platform to solve Linear Programming problems using Simplex, Dual Simplex, and Two-Phase Methods — with real-time visualizations and sensitivity analysis.</i></p>



<br><br>

<a href="#features">✨ Features</a> •
<a href="#demo">📸 Demo</a> •
<a href="#algorithms">📘 Algorithms</a> •
<a href="#installation">⚙️ Installation</a> •
<a href="#usage">🧠 Usage</a>

</div>

<hr>

<h2 id="features">🚀 Key Features</h2>

<ul>
  <li><b>Simplex Method</b> — step-by-step pivots and dictionary updates</li>
  <li><b>Dual Simplex Method</b> — solve infeasible solutions quickly</li>
  <li><b>Two-Phase Simplex</b> — handles ≥, = constraints, and artificial variables</li>
  <li><b>2D Feasible Region Visualization</b> (when n = 2)</li>
  <li><b>Sensitivity Analysis</b> — shadow prices, RHS ranges, allowable changes</li>
  <li><b>Optimal Dictionary View</b> for final basis and reduced costs</li>
  <li><b>Iteration History</b> — tableau sequence for each pivot</li>
  <li><b>Example Problem Library</b> — instantly load standard LP models</li>
  <li><b>Import / Export</b> LP data in JSON</li>
  <li><b>Dark Mode UI</b> with a professional analytics dashboard theme</li>
</ul>

<hr>

<h2 id="demo">📸 Demo Preview</h2>

<p align="center">
  <img width="1004" height="625" alt="Screenshot 2025-11-30 at 8 04 14 PM" src="https://github.com/user-attachments/assets/454fccc5-06a7-4e75-a9f2-1eff0fdff129" />
  <img width="473" height="583" alt="Screenshot 2025-11-30 at 8 04 26 PM" src="https://github.com/user-attachments/assets/62265544-c2ea-44f3-92d2-de6e24656c69" />
  <img width="513" height="636" alt="Screenshot 2025-11-30 at 8 04 36 PM" src="https://github.com/user-attachments/assets/caf5ca57-13c6-41b7-aa02-528b1ed48be6" />



</p>

<hr>

<h2 id="algorithms">📘 Supported Algorithms</h2>

<h3>✔ Simplex Method</h3>
<p>Solves LPs of the form <i>maximize cᵀx subject to Ax ≤ b, x ≥ 0</i>. Shows pivot steps, reduced costs, basis transitions.</p>

<h3>✔ Dual Simplex Method</h3>
<p>Ideal when the initial solution violates feasibility. Automatically selects entering/leaving variables based on dual feasibility.</p>

<h3>✔ Two-Phase Simplex (Phase I + Phase II)</h3>
<ul>
  <li>Handles equality constraints</li>
  <li>Handles ≥ constraints with surplus variables</li>
  <li>Automatically introduces artificial variables</li>
  <li>Phase I finds a BFS; Phase II maximizes/minimizes</li>
</ul>

<hr>

<h2 id="installation">⚙️ Installation</h2>

<pre>
git clone https://github.com/wasayfaizan/linear-optimization-studio
cd linear-optimization-studio
</pre>

Open <code>index.html</code> in your browser — that’s it.

<hr>

<h2 id="usage">🧠 How to Use</h2>

<ol>
  <li>Enter number of constraints (m) and variables (n)</li>
  <li>Fill matrices A, b, c</li>
  <li>Select objective type (Maximize / Minimize)</li>
  <li>Choose constraint types (≤, ≥, =)</li>
  <li>Click <b>Solve LP</b></li>
  <li>Explore:
    <ul>
      <li>Optimal solution</li>
      <li>Slack variables</li>
      <li>Shadow prices</li>
      <li>Optimal dictionary</li>
      <li>Feasible region plot</li>
    </ul>
  </li>
</ol>

<hr>

<h2>🧮 Example Problem</h2>

<pre>
Maximize   60x₁ + 40x₂
Subject to:
  x₁ + x₂ ≤ 6
  2x₁ + x₂ ≤ 9
  2x₁ + 3x₂ ≤ 16
  x₁, x₂ ≥ 0
</pre>

<hr>

<h2>📄 License</h2>
<p>MIT License — free to use, modify, and distribute.</p>

<hr>

<div align="center">
  <h3>Made with ❤️ by Abdul Wasay Faizan</h3>
</div>
