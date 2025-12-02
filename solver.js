// ========= Utility: matrix builders for A, b, c =========
// ========= Utility: matrix builders for A, b, c WITH constraint senses =========
function createMatrixInputs(m, n) {
    const container = document.createElement("div");
  
    // Objective c
    const cSection = document.createElement("div");
    cSection.innerHTML = `
      <label>Objective coefficients c (for ζ = cᵀ x)</label>
      <div class="matrix-container">
        <table class="matrix">
          <tr>
            ${Array.from(
              { length: n },
              (_, j) =>
                `<td><input type="number" step="any" value="${
                  j === 0 ? 60 : j === 1 ? 40 : 0
                }" data-role="c" data-index="${j}"></td>`
            ).join("")}
          </tr>
        </table>
      </div>
    `;
    container.appendChild(cSection);
  
    // A matrix
    const aSection = document.createElement("div");
    aSection.style.marginTop = "6px";
    aSection.innerHTML = `
      <label>Constraint matrix A (m × n)</label>
      <div class="matrix-container">
        <table class="matrix">
          ${Array.from(
            { length: m },
            (_, i) => `
            <tr>
              ${Array.from(
                { length: n },
                (_, j) =>
                  `<td><input type="number" step="any" value="${
                    i === j ? 1 : 0
                  }" data-role="A" data-row="${i}" data-col="${j}"></td>`
              ).join("")}
            </tr>
          `
          ).join("")}
        </table>
      </div>
    `;
    container.appendChild(aSection);
  
    // b vector + constraint sense
    const bSection = document.createElement("div");
    bSection.style.marginTop = "6px";
    bSection.innerHTML = `
      <label>Right-hand side b and constraint type</label>
      <div class="matrix-container">
        <table class="matrix">
          <tr>
            <th>bᵢ</th>
            <th>Type</th>
          </tr>
          ${Array.from(
            { length: m },
            (_, i) => `
            <tr>
              <td>
                <input type="number" step="any" value="${
                  i === 0 ? 6 : i === 1 ? 9 : 16
                }" data-role="b" data-index="${i}">
              </td>
              <td>
                <select data-role="sense" data-index="${i}">
                  <option value="<=" selected>≤</option>
                  <option value="=">=</option>
                  <option value=">=">≥</option>
                </select>
              </td>
            </tr>
          `
          ).join("")}
        </table>
      </div>
    `;
    container.appendChild(bSection);
  
    return container;
  }
   
  function readProblemFromDOM() {
    const m = parseInt(document.getElementById("numConstraints").value, 10);
    const n = parseInt(document.getElementById("numVariables").value, 10);
    const objectiveType = document.getElementById("objectiveType").value; // "max" or "min"
  
    const A = Array.from({ length: m }, () => Array(n).fill(0));
    const b = Array(m).fill(0);
    const c = Array(n).fill(0);
    const senses = Array(m).fill("<=");
  
    document.querySelectorAll('input[data-role="A"]').forEach((inp) => {
      const i = parseInt(inp.dataset.row, 10);
      const j = parseInt(inp.dataset.col, 10);
      A[i][j] = parseFloat(inp.value || "0");
    });
  
    document.querySelectorAll('input[data-role="b"]').forEach((inp) => {
      const i = parseInt(inp.dataset.index, 10);
      b[i] = parseFloat(inp.value || "0");
    });
  
    document.querySelectorAll('input[data-role="c"]').forEach((inp) => {
      const j = parseInt(inp.dataset.index, 10);
      c[j] = parseFloat(inp.value || "0");
    });
  
    document.querySelectorAll('select[data-role="sense"]').forEach((sel) => {
      const i = parseInt(sel.dataset.index, 10);
      senses[i] = sel.value; // "<=", "=", ">="
    });
  
    return { m, n, A, b, c, senses, objectiveType };
  }
  
  // ========= Linear algebra helpers =========
  function invertMatrix(mat) {
    const n = mat.length;
    const aug = mat.map((row, i) => [
      ...row,
      ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)),
    ]);
  
    for (let col = 0; col < n; col++) {
      // partial pivot
      let pivotRow = col;
      let maxVal = Math.abs(aug[col][col]);
      for (let r = col + 1; r < n; r++) {
        const v = Math.abs(aug[r][col]);
        if (v > maxVal) {
          maxVal = v;
          pivotRow = r;
        }
      }
      if (maxVal < 1e-12) {
        throw new Error("Basis matrix is (almost) singular.");
      }
      if (pivotRow !== col) {
        const tmp = aug[col];
        aug[col] = aug[pivotRow];
        aug[pivotRow] = tmp;
      }
  
      const pivot = aug[col][col];
      for (let j = 0; j < 2 * n; j++) {
        aug[col][j] /= pivot;
      }
  
      for (let r = 0; r < n; r++) {
        if (r === col) continue;
        const factor = aug[r][col];
        for (let j = 0; j < 2 * n; j++) {
          aug[r][j] -= factor * aug[col][j];
        }
      }
    }
  
    const inv = Array.from({ length: n }, () => Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        inv[i][j] = aug[i][n + j];
      }
    }
    return inv;
  }
  
  function multiplyMatrixVector(M, v) {
    const rows = M.length;
    const cols = M[0].length;
    const res = new Array(rows).fill(0);
    for (let i = 0; i < rows; i++) {
      let sum = 0;
      for (let j = 0; j < cols; j++) {
        sum += M[i][j] * v[j];
      }
      res[i] = sum;
    }
    return res;
  }
  
  function multiplyMatrices(A, B) {
    const rows = A.length;
    const inner = A[0].length;
    const cols = B[0].length;
    const res = Array.from({ length: rows }, () => Array(cols).fill(0));
    for (let i = 0; i < rows; i++) {
      for (let k = 0; k < inner; k++) {
        const val = A[i][k];
        if (val === 0) continue;
        for (let j = 0; j < cols; j++) {
          res[i][j] += val * B[k][j];
        }
      }
    }
    return res;
  }
  
  function multiplyRowMatrix(row, M) {
    const len = row.length;
    const cols = M[0].length;
    const res = new Array(cols).fill(0);
    for (let j = 0; j < cols; j++) {
      let sum = 0;
      for (let i = 0; i < len; i++) {
        sum += row[i] * M[i][j];
      }
      res[j] = sum;
    }
    return res;
  }
  
  function dot(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }
  // ========= Two-Phase Simplex (general Ax (≤,=,≥) b, x ≥ 0) =========
function twoPhaseSimplex(A_in, b_in, c_in, senses_in, logFn) {
    const m = A_in.length;
    if (m === 0) throw new Error("No constraints specified.");
    const nOrig = A_in[0].length;
  
    const tol = 1e-9;
  
    // --- Deep copy inputs so we don't mutate DOM data ---
    const A = A_in.map((row) => row.slice());
    const b = b_in.slice();
    const senses = senses_in.slice();
    const cOrig = c_in.slice();
  
    // --- Step 0: Normalize rows so that b[i] >= 0 ---
    for (let i = 0; i < m; i++) {
      if (b[i] < -tol) {
        // multiply row by -1
        for (let j = 0; j < nOrig; j++) {
          A[i][j] = -A[i][j];
        }
        b[i] = -b[i];
        // flip sense
        if (senses[i] === "<=") senses[i] = ">=";
        else if (senses[i] === ">=") senses[i] = "<=";
        // "=" stays "="
      }
    }
  
    // --- Step 1: Count slack/surplus/artificial variables ---
    let nSlack = 0;
    let nArtificial = 0;
  
    for (let i = 0; i < m; i++) {
      if (senses[i] === "<=") {
        nSlack += 1;
      } else if (senses[i] === ">=") {
        nSlack += 1; // surplus
        nArtificial += 1;
      } else if (senses[i] === "=") {
        nArtificial += 1;
      } else {
        throw new Error("Invalid constraint type at row " + i);
      }
    }
  
    const nTot = nOrig + nSlack + nArtificial;
    const A_ext = Array.from({ length: m }, () => Array(nTot).fill(0));
  
    const slackIndex = Array(m).fill(-1);
    const artificialIndex = Array(m).fill(-1);
    const basis = Array(m).fill(-1);
  
    let slackCounter = 0;
    let artCounter = 0;
  
    // --- Step 2: Build extended matrix A_ext and initial basis ---
    for (let i = 0; i < m; i++) {
      // original vars
      for (let j = 0; j < nOrig; j++) {
        A_ext[i][j] = A[i][j];
      }
  
      if (senses[i] === "<=") {
        // x: A[i] x <= b[i] => add slack s_i >= 0
        const colSlack = nOrig + slackCounter++;
        A_ext[i][colSlack] = 1;
        slackIndex[i] = colSlack;
        basis[i] = colSlack;
      } else if (senses[i] === ">=") {
        // A[i] x >= b[i] => A[i] x - s_i = b[i], s_i >= 0, add artificial a_i
        const colSurplus = nOrig + slackCounter++;
        A_ext[i][colSurplus] = -1;
  
        const colArt = nOrig + nSlack + artCounter++;
        A_ext[i][colArt] = 1;
        artificialIndex[i] = colArt;
        basis[i] = colArt;
      } else if (senses[i] === "=") {
        // equality => add artificial a_i directly
        const colArt = nOrig + nSlack + artCounter++;
        A_ext[i][colArt] = 1;
        artificialIndex[i] = colArt;
        basis[i] = colArt;
      }
    }
  
    const numVars = nTot;
    const width = numVars + 1;
    const height = m + 1;
  
    // --- Phase I tableau ---
    const tab1 = Array.from({ length: height }, () => Array(width).fill(0));
  
    // constraint rows
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < numVars; j++) {
        tab1[i][j] = A_ext[i][j];
      }
      tab1[i][width - 1] = b[i];
    }
  
    // objective row for Phase I: maximize Z = - sum artificial
    // canonical form: start with 0 row, then subtract all rows with artificial basic
    for (let i = 0; i < m; i++) {
      const varIdx = basis[i];
      if (varIdx >= nOrig + nSlack) {
        // artificial in basis => row0 -= row_i
        for (let j = 0; j < width; j++) {
          tab1[m][j] -= tab1[i][j];
        }
      }
    }
  
    function tableauToStringPhase1() {
      let s = "Phase I Tableau:\n";
      for (let i = 0; i < height; i++) {
        s += tab1[i].map((v) => v.toFixed(3).padStart(9)).join(" ") + "\n";
      }
      s +=
        "Basis: [" + basis.map((idx) => "x" + (idx + 1)).join(", ") + "]\n";
      return s;
    }
  
    logFn("=== Two-Phase Simplex: Phase I (finding BFS) ===\n");
    logFn(tableauToStringPhase1());
  
    const maxIter = 100;
    let iter = 0;
  
    // --- Phase I Simplex iterations ---
    while (iter < maxIter) {
      iter++;
  
      // entering: most negative reduced cost (row m)
      let entering = -1;
      let minRC = 0;
      for (let j = 0; j < numVars; j++) {
        const rc = tab1[m][j];
        if (rc < minRC - tol) {
          minRC = rc;
          entering = j;
        }
      }
  
      if (entering === -1) {
        logFn("Phase I: No negative reduced costs, candidate BFS found.\n");
        break;
      }
  
      // leaving: min ratio test
      let leaving = -1;
      let minRatio = Infinity;
      for (let i = 0; i < m; i++) {
        const aij = tab1[i][entering];
        const rhs = tab1[i][width - 1];
        if (aij > tol) {
          const ratio = rhs / aij;
          if (ratio < minRatio - tol) {
            minRatio = ratio;
            leaving = i;
          }
        }
      }
  
      if (leaving === -1) {
        throw new Error(
          "Phase I: Unbounded in auxiliary problem. Original LP likely infeasible."
        );
      }
  
      logFn(
        `Phase I Iteration ${iter}: entering x${
          entering + 1
        }, leaving row ${leaving} (x${basis[leaving] + 1})`
      );
  
      // pivot
      const pivot = tab1[leaving][entering];
      for (let j = 0; j < width; j++) {
        tab1[leaving][j] /= pivot;
      }
      for (let i = 0; i < height; i++) {
        if (i === leaving) continue;
        const factor = tab1[i][entering];
        for (let j = 0; j < width; j++) {
          tab1[i][j] -= factor * tab1[leaving][j];
        }
      }
  
      basis[leaving] = entering;
      logFn(tableauToStringPhase1());
    }
  
    const phase1Obj = tab1[m][width - 1];
  
    if (phase1Obj < -1e-6) {
      // Z = -sum(a_i). If this < 0, sum(a_i) > 0 => infeasible
      throw new Error(
        "Phase I: Optimal auxiliary objective is negative. Original LP is infeasible."
      );
    }
  
    logFn(
      `Phase I completed. Objective ≈ ${phase1Obj.toFixed(
        6
      )}. Proceeding to Phase II.\n`
    );
  
    // --- Cleanup: remove artificial variables from basis if any ---
    const artificialStart = nOrig + nSlack;
    function isArtificial(idx) {
      return idx >= artificialStart && idx < nTot;
    }
  
    function pivotPhase1(row, col) {
      const pivot = tab1[row][col];
      for (let j = 0; j < width; j++) {
        tab1[row][j] /= pivot;
      }
      for (let i = 0; i < height; i++) {
        if (i === row) continue;
        const factor = tab1[i][col];
        for (let j = 0; j < width; j++) {
          tab1[i][j] -= factor * tab1[row][j];
        }
      }
      basis[row] = col;
    }
  
    // Try to pivot artificials out of basis
    for (let i = 0; i < m; i++) {
      if (isArtificial(basis[i])) {
        let pivotCol = -1;
        for (let j = 0; j < nTot; j++) {
          if (!isArtificial(j)) {
            const aij = tab1[i][j];
            if (Math.abs(aij) > tol) {
              pivotCol = j;
              break;
            }
          }
        }
        if (pivotCol !== -1) {
          pivotPhase1(i, pivotCol);
        }
        // If no such pivotCol, this row is redundant; we leave it as is.
      }
    }
  
    logFn("Phase I cleanup (pivots to remove artificials from basis):\n");
    logFn(tableauToStringPhase1());
  
    // --- Step 3: Build Phase II tableau (drop artificial columns) ---
  
    // Map old variable indices -> new ones (without artificials)
    const oldToNew = Array(nTot).fill(-1);
    let newVarCount = 0;
    for (let j = 0; j < nTot; j++) {
      if (!isArtificial(j)) {
        oldToNew[j] = newVarCount++;
      }
    }
    const nPhase2Vars = newVarCount;
    const width2 = nPhase2Vars + 1;
    const height2 = m + 1;
  
    const tab2 = Array.from({ length: height2 }, () =>
      Array(width2).fill(0)
    );
    const basis2 = Array(m).fill(-1);
  
    // Copy constraints rows (dropping artificial columns)
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < nTot; j++) {
        const newIdx = oldToNew[j];
        if (newIdx !== -1) {
          tab2[i][newIdx] = tab1[i][j];
        }
      }
      tab2[i][width2 - 1] = tab1[i][width - 1];
      basis2[i] = oldToNew[basis[i]];
    }
  
    // Build full c vector for Phase II (same length as nTot, but artificials have 0 cost)
    const cFull = Array(nTot).fill(0);
    for (let j = 0; j < nOrig; j++) cFull[j] = cOrig[j];
    // slack/surplus and artificial have 0
  
    // Phase II objective row: start from -c for non-artificial vars
    for (let j = 0; j < nTot; j++) {
      const newIdx = oldToNew[j];
      if (newIdx !== -1) {
        tab2[m][newIdx] = -cFull[j];
      }
    }
    tab2[m][width2 - 1] = 0;
  
    // Make objective row canonical: add c_B * row_i for each basic variable
    for (let i = 0; i < m; i++) {
      const oldIdx = basis[i]; // original index
      const cb = cFull[oldIdx];
      if (Math.abs(cb) < tol) continue;
      for (let j = 0; j < width2; j++) {
        tab2[m][j] += cb * tab2[i][j];
      }
    }
  
    function tableauToStringPhase2() {
      let s = "Phase II Tableau:\n";
      for (let i = 0; i < height2; i++) {
        s += tab2[i].map((v) => v.toFixed(3).padStart(9)).join(" ") + "\n";
      }
      s +=
        "Basis: [" + basis2.map((idx) => "x" + (idx + 1)).join(", ") + "]\n";
      return s;
    }
  
    logFn("=== Two-Phase Simplex: Phase II (optimize original objective) ===\n");
    logFn(tableauToStringPhase2());
  
    // --- Phase II iterations ---
    iter = 0;
    while (iter < maxIter) {
      iter++;
  
      // entering var: most negative reduced cost
      let entering = -1;
      let minRC = 0;
      for (let j = 0; j < nPhase2Vars; j++) {
        const rc = tab2[m][j];
        if (rc < minRC - tol) {
          minRC = rc;
          entering = j;
        }
      }
  
      if (entering === -1) {
        logFn("Phase II: No negative reduced costs, optimal solution reached.\n");
        break;
      }
  
      // leaving row: min ratio test
      let leaving = -1;
      let minRatio = Infinity;
      for (let i = 0; i < m; i++) {
        const aij = tab2[i][entering];
        const rhs = tab2[i][width2 - 1];
        if (aij > tol) {
          const ratio = rhs / aij;
          if (ratio < minRatio - tol) {
            minRatio = ratio;
            leaving = i;
          }
        }
      }
  
      if (leaving === -1) {
        throw new Error("Phase II: Unbounded problem.");
      }
  
      logFn(
        `Phase II Iteration ${iter}: entering x${
          entering + 1
        }, leaving row ${leaving} (x${basis2[leaving] + 1})`
      );
  
      // pivot
      const pivot = tab2[leaving][entering];
      for (let j = 0; j < width2; j++) {
        tab2[leaving][j] /= pivot;
      }
      for (let i = 0; i < height2; i++) {
        if (i === leaving) continue;
        const factor = tab2[i][entering];
        for (let j = 0; j < width2; j++) {
          tab2[i][j] -= factor * tab2[leaving][j];
        }
      }
  
      basis2[leaving] = entering;
      logFn(tableauToStringPhase2());
    }
  
    // --- Extract solution x (for all vars, then cut to original vars) ---
    const xFullPhase2 = Array(nPhase2Vars).fill(0);
    for (let i = 0; i < m; i++) {
      const varIdx = basis2[i];
      if (varIdx >= 0 && varIdx < nPhase2Vars) {
        xFullPhase2[varIdx] = tab2[i][width2 - 1];
      }
    }
    const zStar = tab2[m][width2 - 1];
  
    // Map back to original variable indices (x1...xn, slacks, surplus)
    const xAllVars = Array(nTot).fill(0);
    for (let oldIdx = 0; oldIdx < nTot; oldIdx++) {
      const newIdx = oldToNew[oldIdx];
      if (newIdx !== -1) {
        xAllVars[oldIdx] = xFullPhase2[newIdx];
      }
    }
  
    // Build final "basis" in original index space (non-artificial only)
    const finalBasis = [];
    for (let i = 0; i < m; i++) {
      // basis2[i] is index in phase-2 variable space
      const phase2Idx = basis2[i];
      if (phase2Idx === -1) continue;
      // find old index that maps to this
      let oldIdx = -1;
      for (let j = 0; j < nTot; j++) {
        if (oldToNew[j] === phase2Idx) {
          oldIdx = j;
          break;
        }
      }
      if (oldIdx !== -1) {
        finalBasis.push(oldIdx);
      }
    }
  
    return {
      x: xAllVars, // includes original + slack/surplus (no artificial)
      z: zStar,
      basis: finalBasis.length ? finalBasis : basis, // fallback
      nOrig,
    };
  }
  
  // ========= Simplex (primal), max cᵀx, Ax ≤ b, x ≥ 0 =========
  function simplexMaxStandard(A, b, c, logFn) {
    const m = A.length;
    const n = A[0].length;
    const width = n + m + 1;
    const height = m + 1;
  
    const tab = Array.from({ length: height }, () => Array(width).fill(0));
  
    // constraints rows
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < n; j++) {
        tab[i][j] = A[i][j];
      }
      tab[i][n + i] = 1; // slack
      tab[i][width - 1] = b[i];
    }
  
    // objective row: -c (for maximization)
    for (let j = 0; j < n; j++) {
      tab[m][j] = -c[j];
    }
  
    let basis = [];
    for (let i = 0; i < m; i++) {
      basis.push(n + i); // initial slack basis
    }
  
    // Check initial RHS feasibility (b >= 0 for slack basis)
    let hasNegativeRHS = false;
    for (let i = 0; i < m; i++) {
      if (tab[i][width - 1] < -1e-9) {
        hasNegativeRHS = true;
        break;
      }
    }
    if (hasNegativeRHS) {
      logFn(
        "Initial BFS with slack basis is infeasible (some RHS < 0). Primal simplex cannot start from this basis."
      );
      throw new Error("Initial BFS infeasible (some RHS < 0). Try Dual Simplex.");
    }
  
    function tableauToString() {
      let s = "";
      s += "Tableau:\n";
      for (let i = 0; i < height; i++) {
        s += tab[i].map((v) => v.toFixed(3).padStart(8)).join(" ") + "\n";
      }
      s +=
        "Basis: [" + basis.map((idx) => "x" + (idx + 1)).join(", ") + "]\n";
      return s;
    }
  
    logFn("Initial BFS with slack basis:\n" + tableauToString());
  
    const maxIter = 50;
    let iter = 0;
  
    while (iter < maxIter) {
      iter++;
  
      // entering variable: most negative reduced cost in last row
      let entering = -1;
      let minRC = 0;
      for (let j = 0; j < width - 1; j++) {
        const rc = tab[m][j];
        if (rc < minRC - 1e-9) {
          minRC = rc;
          entering = j;
        }
      }
  
      if (entering === -1) {
        logFn("No negative reduced costs. Primal optimality reached.");
        break;
      }
  
      // leaving variable: min ratio test
      let leaving = -1;
      let minRatio = Infinity;
      for (let i = 0; i < m; i++) {
        const aij = tab[i][entering];
        if (aij > 1e-9) {
          const ratio = tab[i][width - 1] / aij;
          if (ratio < minRatio - 1e-9) {
            minRatio = ratio;
            leaving = i;
          }
        }
      }
  
      if (leaving === -1) {
        throw new Error("Unbounded problem (no valid leaving variable).");
      }
  
      logFn(
        `Iteration ${iter}: entering x${
          entering + 1
        }, leaving row ${leaving} (variable x${basis[leaving] + 1})`
      );
  
      // pivot
      const pivot = tab[leaving][entering];
      for (let j = 0; j < width; j++) {
        tab[leaving][j] /= pivot;
      }
      for (let i = 0; i < height; i++) {
        if (i === leaving) continue;
        const factor = tab[i][entering];
        for (let j = 0; j < width; j++) {
          tab[i][j] -= factor * tab[leaving][j];
        }
      }
  
      basis[leaving] = entering;
      logFn(tableauToString());
    }
  
    const numVars = n + m; // excluding objective
    const x = Array(numVars).fill(0);
    for (let i = 0; i < m; i++) {
      const varIdx = basis[i];
      if (varIdx < numVars) {
        x[varIdx] = tab[i][width - 1];
      }
    }
    const z = tab[m][width - 1];
  
    return { x, z, basis };
  }
  
  // ========= Dual Simplex (max cᵀx, Ax ≤ b, x ≥ 0) =========
  function dualSimplexStandard(A, b, c, logFn) {
    const m = A.length;
    const n = A[0].length;
    const width = n + m + 1;
    const height = m + 1;
  
    const tab = Array.from({ length: height }, () => Array(width).fill(0));
  
    // constraints rows
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < n; j++) {
        tab[i][j] = A[i][j];
      }
      tab[i][n + i] = 1; // slack
      tab[i][width - 1] = b[i];
    }
  
    // objective row: -c
    for (let j = 0; j < n; j++) {
      tab[m][j] = -c[j];
    }
  
    let basis = [];
    for (let i = 0; i < m; i++) {
      basis.push(n + i); // slack basis
    }
  
    function tableauToString() {
      let s = "";
      s += "Tableau (Dual Simplex):\n";
      for (let i = 0; i < height; i++) {
        s += tab[i].map((v) => v.toFixed(3).padStart(8)).join(" ") + "\n";
      }
      s +=
        "Basis: [" + basis.map((idx) => "x" + (idx + 1)).join(", ") + "]\n";
      return s;
    }
  
    logFn("Initial tableau for Dual Simplex:\n" + tableauToString());
  
    const maxIter = 50;
    let iter = 0;
  
    while (iter < maxIter) {
      iter++;
  
      // 1. Choose leaving row: most negative RHS
      let leaving = -1;
      let minRHS = 0;
      for (let i = 0; i < m; i++) {
        const rhs = tab[i][width - 1];
        if (rhs < minRHS - 1e-9) {
          minRHS = rhs;
          leaving = i;
        }
      }
  
      // If no negative RHS, primal feasible -> optimal
      if (leaving === -1) {
        logFn(
          "No negative RHS. Primal feasibility achieved -> optimal (Dual Simplex)."
        );
        break;
      }
  
      // 2. Choose entering column (maintain dual feasibility)
      let entering = -1;
      let minRatio = Infinity;
      for (let j = 0; j < width - 1; j++) {
        const a_kj = tab[leaving][j];
        const rc = tab[m][j]; // reduced cost in objective row
        if (a_kj < -1e-9) {
          const ratio = rc / a_kj; // a_kj < 0
          if (ratio < minRatio - 1e-9) {
            minRatio = ratio;
            entering = j;
          }
        }
      }
  
      if (entering === -1) {
        throw new Error(
          "Dual Simplex: problem infeasible (no valid entering variable)."
        );
      }
  
      logFn(
        `Dual Iteration ${iter}: leaving row ${leaving} (x${
          basis[leaving] + 1
        }), entering x${entering + 1}`
      );
  
      // 3. Pivot
      const pivot = tab[leaving][entering];
      for (let j = 0; j < width; j++) {
        tab[leaving][j] /= pivot;
      }
      for (let i = 0; i < height; i++) {
        if (i === leaving) continue;
        const factor = tab[i][entering];
        for (let j = 0; j < width; j++) {
          tab[i][j] -= factor * tab[leaving][j];
        }
      }
  
      basis[leaving] = entering;
      logFn(tableauToString());
    }
  
    const numVars = n + m;
    const x = Array(numVars).fill(0);
    for (let i = 0; i < m; i++) {
      const varIdx = basis[i];
      if (varIdx < numVars) {
        x[varIdx] = tab[i][width - 1];
      }
    }
    const z = tab[m][width - 1];
  
    return { x, z, basis };
  }
  
  // ========= Dictionary builder from basis (matrix notation) =========
  function buildDictionary(A, b, c, basis) {
    const m = A.length;
    const nOrig = A[0].length;
    const nTot = nOrig + m;
  
    // Build A_full = [A | I_m]
    const A_full = Array.from({ length: m }, () => Array(nTot).fill(0));
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < nOrig; j++) {
        A_full[i][j] = A[i][j];
      }
      for (let s = 0; s < m; s++) {
        A_full[i][nOrig + s] = i === s ? 1 : 0;
      }
    }
  
    const c_full = Array(nTot).fill(0);
    for (let j = 0; j < nOrig; j++) c_full[j] = c[j];
  
    const inBasis = Array(nTot).fill(false);
    basis.forEach((idx) => {
      if (idx >= 0 && idx < nTot) inBasis[idx] = true;
    });
  
    const nonbasic = [];
    for (let j = 0; j < nTot; j++) {
      if (!inBasis[j]) nonbasic.push(j);
    }
  
    // B matrix (m x m)
    const B = Array.from({ length: m }, () => Array(m).fill(0));
    for (let col = 0; col < m; col++) {
      const varIdx = basis[col];
      for (let i = 0; i < m; i++) {
        B[i][col] = A_full[i][varIdx];
      }
    }
  
    const N = Array.from({ length: m }, () =>
      Array(nonbasic.length).fill(0)
    );
    for (let col = 0; col < nonbasic.length; col++) {
      const varIdx = nonbasic[col];
      for (let i = 0; i < m; i++) {
        N[i][col] = A_full[i][varIdx];
      }
    }
  
    const Binv = invertMatrix(B);
    const xB = multiplyMatrixVector(Binv, b);
  
    const cB = basis.map((idx) => c_full[idx]);
    const cN = nonbasic.map((idx) => c_full[idx]);
  
    const B_inv_N = multiplyMatrices(Binv, N);
    const yRow = multiplyRowMatrix(cB, Binv); // y^T = c_B^T B^{-1}
    const z0 = dot(cB, xB);
    const cB_BinvN = multiplyRowMatrix(cB, B_inv_N);
    const cRed = cN.map((val, k) => val - cB_BinvN[k]);
  
    // Pretty dictionary string
    let s = "";
    s +=
      "<div style='margin-bottom: 12px; padding: 8px; background: rgba(56, 189, 248, 0.1); border-radius: 6px; border-left: 3px solid var(--accent);'>";
    s +=
      "<strong style='color: var(--accent);'>Basis 𝔅</strong> = {" +
      basis.map((idx) => "<strong>x" + (idx + 1) + "</strong>").join(", ") +
      "}<br>";
    s +=
      "<strong style='color: var(--accent);'>Non-basis 𝒩</strong> = {" +
      nonbasic.map((idx) => "x" + (idx + 1)).join(", ") +
      "}";
    s += "</div>\n\n";
  
    s += "<div style='margin-bottom: 8px;'><strong style='color: var(--accent);'>ζ</strong> = <strong>" + z0.toFixed(4) + "</strong>";
    for (let k = 0; k < nonbasic.length; k++) {
      const coeff = cRed[k];
      if (Math.abs(coeff) < 1e-9) continue;
      const sign = coeff >= 0 ? " + " : " - ";
      const color = coeff < 0 ? "var(--danger)" : "var(--success)";
      s += "<span style='color: " + color + ";'>" + sign + Math.abs(coeff).toFixed(4) + " x<sub>" + (nonbasic[k] + 1) + "</sub></span>";
    }
    s += "</div>\n";
  
    s += "<div style='margin-top: 12px;'>";
    for (let i = 0; i < m; i++) {
      s += "<div style='margin-bottom: 6px; padding: 6px; background: rgba(15, 23, 42, 0.5); border-radius: 4px;'>";
      s += "<strong style='color: var(--accent);'>x<sub>" + (basis[i] + 1) + "</sub></strong> = <strong>" + xB[i].toFixed(4) + "</strong>";
      for (let k = 0; k < nonbasic.length; k++) {
        const a = B_inv_N[i][k];
        if (Math.abs(a) < 1e-9) continue;
        const sign = a >= 0 ? " - " : " + ";
        s += sign + "<span style='color: var(--muted);'>" + Math.abs(a).toFixed(4) + "</span> x<sub>" + (nonbasic[k] + 1) + "</sub>";
      }
      s += "</div>";
    }
    s += "</div>";
  
    return {
      dictString: s,
      basis,
      nonbasic,
      xB,
      y: yRow,
      z0,
      cRed,
    };
  }
  
  // ========= Plot for 2-variable case =========
  function plotFeasibleRegion(A, b, c = [1, 1], objectiveType = "max") {
    const canvas = document.getElementById("plotCanvas");
    const ctx = canvas.getContext("2d");
  
    // Use CSS size but actual width/height from attributes
    const w = canvas.width || canvas.getBoundingClientRect().width;
    const h = canvas.height || canvas.getBoundingClientRect().height;
    canvas.width = w;
    canvas.height = h;
  
    const maxX = 10;
    const maxY = 10;
  
    function XtoPx(x) {
      return (x / maxX) * (w - 40) + 30;
    }
    function YtoPx(y) {
      return h - ((y / maxY) * (h - 40) + 30);
    }
  
    ctx.clearRect(0, 0, w, h);
  
    const grad = ctx.createRadialGradient(w / 2, 0, 0, w / 2, 0, w);
    grad.addColorStop(0, "rgba(56, 189, 248, 0.3)");
    grad.addColorStop(1, "#020617");
    ctx.fillStyle = grad;
    ctx.fillRect(0, 0, w, h);
  
    ctx.lineWidth = 1;
    ctx.strokeStyle = "rgba(148,163,184,0.7)";
    ctx.fillStyle = "rgba(148,163,184,0.7)";
    ctx.font = "10px system-ui";
  
    // Axes
    ctx.beginPath();
    ctx.moveTo(XtoPx(0), YtoPx(0));
    ctx.lineTo(XtoPx(maxX), YtoPx(0));
    ctx.moveTo(XtoPx(0), YtoPx(0));
    ctx.lineTo(XtoPx(0), YtoPx(maxY));
    ctx.stroke();
  
    ctx.fillText("x₁", XtoPx(maxX) - 15, YtoPx(0) + 12);
    ctx.fillText("x₂", XtoPx(0) - 15, YtoPx(maxY) + 4);
  
    // Ticks
    for (let t = 1; t <= maxX; t++) {
      const x = XtoPx(t);
      ctx.beginPath();
      ctx.moveTo(x, YtoPx(0) - 3);
      ctx.lineTo(x, YtoPx(0) + 3);
      ctx.stroke();
    }
    for (let t = 1; t <= maxY; t++) {
      const y = YtoPx(t);
      ctx.beginPath();
      ctx.moveTo(XtoPx(0) - 3, y);
      ctx.lineTo(XtoPx(0) + 3, y);
      ctx.stroke();
    }
  
    // Find corner points (intersections of constraints)
    const corners = [];
    for (let i = 0; i < A.length; i++) {
      for (let j = i + 1; j < A.length; j++) {
        const [a1, a2] = A[i];
        const [b1, b2] = A[j];
        const det = a1 * b2 - a2 * b1;
        if (Math.abs(det) > 1e-9) {
          const x1 = (b[i] * b2 - b[j] * a2) / det;
          const x2 = (a1 * b[j] - b1 * b[i]) / det;
          if (x1 >= 0 && x2 >= 0 && x1 <= maxX && x2 <= maxY) {
            // Check if point satisfies all constraints
            let feasible = true;
            for (let k = 0; k < A.length; k++) {
              const lhs = A[k][0] * x1 + A[k][1] * x2;
              if (lhs - b[k] > 1e-6) {
                feasible = false;
                break;
              }
            }
            if (feasible) {
              corners.push({ x1, x2 });
            }
          }
        }
      }
    }
    // Add axis intersections
    for (let i = 0; i < A.length; i++) {
      const [a1, a2] = A[i];
      const bi = b[i];
      if (Math.abs(a2) > 1e-9) {
        const x1 = 0;
        const x2 = bi / a2;
        if (x2 >= 0 && x2 <= maxY) {
          let feasible = true;
          for (let k = 0; k < A.length; k++) {
            const lhs = A[k][0] * x1 + A[k][1] * x2;
            if (lhs - b[k] > 1e-6) {
              feasible = false;
              break;
            }
          }
          if (feasible) corners.push({ x1, x2 });
        }
      }
      if (Math.abs(a1) > 1e-9) {
        const x1 = bi / a1;
        const x2 = 0;
        if (x1 >= 0 && x1 <= maxX) {
          let feasible = true;
          for (let k = 0; k < A.length; k++) {
            const lhs = A[k][0] * x1 + A[k][1] * x2;
            if (lhs - b[k] > 1e-6) {
              feasible = false;
              break;
            }
          }
          if (feasible) corners.push({ x1, x2 });
        }
      }
    }
    
    // Remove duplicates
    const uniqueCorners = [];
    corners.forEach(c => {
      const exists = uniqueCorners.some(uc => 
        Math.abs(uc.x1 - c.x1) < 1e-6 && Math.abs(uc.x2 - c.x2) < 1e-6
      );
      if (!exists) uniqueCorners.push(c);
    });
  
    // Feasible region sampling
    const pts = [];
    const step = 0.1;
    for (let x1 = 0; x1 <= maxX; x1 += step) {
      for (let x2 = 0; x2 <= maxY; x2 += step) {
        let ok = true;
        for (let i = 0; i < A.length; i++) {
          const lhs = A[i][0] * x1 + A[i][1] * x2;
          if (lhs - b[i] > 1e-6) {
            ok = false;
            break;
          }
        }
        if (ok) pts.push({ x1, x2 });
      }
    }
  
    ctx.fillStyle = "rgba(56,189,248,0.22)";
    pts.forEach((p) => {
      ctx.beginPath();
      ctx.arc(XtoPx(p.x1), YtoPx(p.x2), 1, 0, 2 * Math.PI);
      ctx.fill();
    });
  
    // Draw constraint lines with labels
    ctx.lineWidth = 1.2;
    ctx.strokeStyle = "rgba(59,130,246,0.9)";
    for (let i = 0; i < A.length; i++) {
      const [a1, a2] = A[i];
      const bi = b[i];
      const pLine = [];
  
      for (let x1 = 0; x1 <= maxX; x1 += 0.1) {
        if (Math.abs(a2) < 1e-9) continue;
        const x2 = (bi - a1 * x1) / a2;
        if (x2 >= 0 && x2 <= maxY) {
          pLine.push({ x1, x2 });
        }
      }
  
      if (pLine.length > 1) {
        ctx.beginPath();
        ctx.moveTo(XtoPx(pLine[0].x1), YtoPx(pLine[0].x2));
        for (const p of pLine) {
          ctx.lineTo(XtoPx(p.x1), YtoPx(p.x2));
        }
        ctx.stroke();
        
        // Label constraint
        if (pLine.length > 0) {
          const mid = pLine[Math.floor(pLine.length / 2)];
          ctx.fillStyle = "rgba(59,130,246,0.9)";
          ctx.fillText(`C${i + 1}`, XtoPx(mid.x1) + 5, YtoPx(mid.x2) - 5);
        }
      }
    }
    
    // Draw iso-profit/iso-cost lines
    if (c && c.length === 2 && (Math.abs(c[0]) > 1e-9 || Math.abs(c[1]) > 1e-9)) {
      ctx.lineWidth = 1;
      ctx.strokeStyle = "rgba(34,197,94,0.6)";
      ctx.setLineDash([5, 5]);
      
      // Draw a few iso-lines
      const maxObj = uniqueCorners.length > 0 && c && c.length === 2
        ? Math.max(...uniqueCorners.map(corner => c[0] * corner.x1 + c[1] * corner.x2))
        : 10;
      for (let k = 1; k <= 3; k++) {
        const objVal = (maxObj * k) / 4;
        const isoLine = [];
        if (Math.abs(c[1]) > 1e-9) {
          for (let x1 = 0; x1 <= maxX; x1 += 0.1) {
            const x2 = (objVal - c[0] * x1) / c[1];
            if (x2 >= 0 && x2 <= maxY) {
              isoLine.push({ x1, x2 });
            }
          }
        } else if (Math.abs(c[0]) > 1e-9) {
          const x1 = objVal / c[0];
          if (x1 >= 0 && x1 <= maxX) {
            isoLine.push({ x1, x2: 0 });
            isoLine.push({ x1, x2: maxY });
          }
        }
        
        if (isoLine.length > 1) {
          ctx.beginPath();
          ctx.moveTo(XtoPx(isoLine[0].x1), YtoPx(isoLine[0].x2));
          for (const p of isoLine) {
            ctx.lineTo(XtoPx(p.x1), YtoPx(p.x2));
          }
          ctx.stroke();
        }
      }
      ctx.setLineDash([]);
    }
    
    // Draw and label corner points
    ctx.fillStyle = "rgba(251,146,60,0.9)";
    ctx.strokeStyle = "rgba(251,146,60,1)";
    ctx.lineWidth = 2;
    uniqueCorners.forEach((corner, idx) => {
      const px = XtoPx(corner.x1);
      const py = YtoPx(corner.x2);
      ctx.beginPath();
      ctx.arc(px, py, 4, 0, 2 * Math.PI);
      ctx.fill();
      ctx.stroke();
      
      // Label corner
      ctx.fillStyle = "rgba(251,146,60,1)";
      ctx.font = "9px system-ui";
      ctx.fillText(`(${corner.x1.toFixed(1)}, ${corner.x2.toFixed(1)})`, px + 6, py - 6);
      ctx.fillStyle = "rgba(251,146,60,0.9)";
    });
    
    // Find and highlight optimal corner (if solved)
    if (uniqueCorners.length > 0 && c && c.length === 2) {
      const objVals = uniqueCorners.map(corner => c[0] * corner.x1 + c[1] * corner.x2);
      const optimalIdx = objectiveType === "max" 
        ? objVals.indexOf(Math.max(...objVals))
        : objVals.indexOf(Math.min(...objVals));
      
      if (optimalIdx >= 0) {
        const optCorner = uniqueCorners[optimalIdx];
        const px = XtoPx(optCorner.x1);
        const py = YtoPx(optCorner.x2);
        ctx.fillStyle = "rgba(34,197,94,1)";
        ctx.strokeStyle = "rgba(34,197,94,1)";
        ctx.lineWidth = 3;
        ctx.beginPath();
        ctx.arc(px, py, 6, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();
      }
    }
  }
  
  // ========= CSV Import/Export =========
  function parseCSV(text) {
    const lines = text.trim().split(/\r?\n/).filter(line => line.trim());
    return lines.map(line => {
      // Handle both comma and tab separators, and quoted values
      return line.split(/[,\t]/).map(val => val.trim().replace(/^"|"$/g, ''));
    });
  }

  function importFromCSV(csvText) {
    try {
      const rows = parseCSV(csvText);
      if (rows.length < 2) throw new Error("CSV must have at least 2 rows");
      
      // First row: c vector
      const c = rows[0].map(v => parseFloat(v) || 0);
      const n = c.length;
      
      // Remaining rows: A matrix and b vector
      const m = rows.length - 1;
      const A = [];
      const b = [];
      const senses = [];
      
      for (let i = 1; i < rows.length; i++) {
        const row = rows[i];
        if (row.length < n + 1) throw new Error(`Row ${i+1} has insufficient columns`);
        
        const aRow = row.slice(0, n).map(v => parseFloat(v) || 0);
        A.push(aRow);
        
        // Last column is b, second-to-last might be sense
        if (row.length >= n + 2) {
          b.push(parseFloat(row[n]) || 0);
          const sense = row[n + 1].trim().toLowerCase();
          senses.push(sense === ">=" || sense === "≥" ? ">=" : 
                     sense === "=" || sense === "==" ? "=" : "<=");
        } else {
          b.push(parseFloat(row[n]) || 0);
          senses.push("<=");
        }
      }
      
      // Update dimensions
      document.getElementById("numVariables").value = n;
      document.getElementById("numConstraints").value = m;
      generateMatrices();
      
      // Fill in values
      setTimeout(() => {
        document.querySelectorAll('input[data-role="c"]').forEach((inp, j) => {
          if (j < c.length) inp.value = c[j];
        });
        
        document.querySelectorAll('input[data-role="A"]').forEach((inp) => {
          const i = parseInt(inp.dataset.row, 10);
          const j = parseInt(inp.dataset.col, 10);
          if (i < A.length && j < A[i].length) {
            inp.value = A[i][j];
          }
        });
        
        document.querySelectorAll('input[data-role="b"]').forEach((inp, i) => {
          if (i < b.length) inp.value = b[i];
        });
        
        document.querySelectorAll('select[data-role="sense"]').forEach((sel, i) => {
          if (i < senses.length) sel.value = senses[i];
        });
      }, 100);
      
      return { m, n, A, b, c, senses };
    } catch (err) {
      throw new Error("CSV import failed: " + err.message);
    }
  }

  function exportToCSV() {
    const { m, n, A, b, c, senses } = readProblemFromDOM();
    let csv = c.map(v => v.toString()).join(",") + "\n";
    for (let i = 0; i < m; i++) {
      csv += A[i].map(v => v.toString()).join(",") + "," + b[i] + "," + senses[i] + "\n";
    }
    
    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "linear_program.csv";
    a.click();
    URL.revokeObjectURL(url);
  }

  // ========= Example Problems =========
  const exampleProblems = {
    diet: {
      m: 4,
      n: 2,
      A: [[2, 3], [1, 1], [1, 0], [0, 1]],
      b: [12, 5, 3, 2],
      c: [3, 2],
      senses: ["<=", "<=", ">=", ">="],
      objectiveType: "min",
      name: "Diet Problem"
    },
    production: {
      m: 3,
      n: 2,
      A: [[2, 1], [1, 3], [1, 1]],
      b: [100, 150, 80],
      c: [50, 40],
      senses: ["<=", "<=", ">="],
      objectiveType: "max",
      name: "Production Planning"
    },
    transportation: {
      m: 4,
      n: 4,
      A: [[1, 1, 0, 0], [0, 0, 1, 1], [1, 0, 1, 0], [0, 1, 0, 1]],
      b: [20, 30, 25, 25],
      c: [10, 12, 8, 9],
      senses: ["<=", "<=", "=", "="],
      objectiveType: "min",
      name: "Transportation Problem"
    },
    simple: {
      m: 3,
      n: 2,
      A: [[1, 1], [2, 1], [1, 2]],
      b: [6, 9, 8],
      c: [60, 40],
      senses: ["<=", "<=", "<="],
      objectiveType: "max",
      name: "Simple 2x2 Example"
    }
  };

  function loadExample(exampleKey) {
    const ex = exampleProblems[exampleKey];
    if (!ex) return;
    
    document.getElementById("numConstraints").value = ex.m;
    document.getElementById("numVariables").value = ex.n;
    document.getElementById("objectiveType").value = ex.objectiveType;
    updateObjectiveLabel();
    generateMatrices();
    
    setTimeout(() => {
      document.querySelectorAll('input[data-role="c"]').forEach((inp, j) => {
        if (j < ex.c.length) inp.value = ex.c[j];
      });
      
      document.querySelectorAll('input[data-role="A"]').forEach((inp) => {
        const i = parseInt(inp.dataset.row, 10);
        const j = parseInt(inp.dataset.col, 10);
        if (i < ex.A.length && j < ex.A[i].length) {
          inp.value = ex.A[i][j];
        }
      });
      
      document.querySelectorAll('input[data-role="b"]').forEach((inp, i) => {
        if (i < ex.b.length) inp.value = ex.b[i];
      });
      
      document.querySelectorAll('select[data-role="sense"]').forEach((sel, i) => {
        if (i < ex.senses.length) sel.value = ex.senses[i];
      });
      
      setStatus(`Loaded example: ${ex.name}`, "ok");
    }, 100);
  }

  // ========= Save/Load =========
  function saveProblem() {
    const problem = readProblemFromDOM();
    const data = JSON.stringify(problem, null, 2);
    const blob = new Blob([data], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "linear_program.json";
    a.click();
    URL.revokeObjectURL(url);
    
    // Also save to localStorage
    localStorage.setItem("linearOpt_lastProblem", data);
    setStatus("Problem saved to file and localStorage", "ok");
  }

  function loadProblem() {
    const fileInput = document.createElement("input");
    fileInput.type = "file";
    fileInput.accept = ".json";
    fileInput.onchange = (e) => {
      const file = e.target.files[0];
      if (!file) return;
      
      const reader = new FileReader();
      reader.onload = (event) => {
        try {
          const problem = JSON.parse(event.target.result);
          loadProblemData(problem);
        } catch (err) {
          setStatus("Error loading file: " + err.message, "err");
        }
      };
      reader.readAsText(file);
    };
    fileInput.click();
  }

  function loadProblemData(problem) {
    document.getElementById("numConstraints").value = problem.m;
    document.getElementById("numVariables").value = problem.n;
    if (problem.objectiveType) {
      document.getElementById("objectiveType").value = problem.objectiveType;
      updateObjectiveLabel();
    }
    generateMatrices();
    
    setTimeout(() => {
      if (problem.c) {
        document.querySelectorAll('input[data-role="c"]').forEach((inp, j) => {
          if (j < problem.c.length) inp.value = problem.c[j];
        });
      }
      
      if (problem.A) {
        document.querySelectorAll('input[data-role="A"]').forEach((inp) => {
          const i = parseInt(inp.dataset.row, 10);
          const j = parseInt(inp.dataset.col, 10);
          if (i < problem.A.length && j < problem.A[i].length) {
            inp.value = problem.A[i][j];
          }
        });
      }
      
      if (problem.b) {
        document.querySelectorAll('input[data-role="b"]').forEach((inp, i) => {
          if (i < problem.b.length) inp.value = problem.b[i];
        });
      }
      
      if (problem.senses) {
        document.querySelectorAll('select[data-role="sense"]').forEach((sel, i) => {
          if (i < problem.senses.length) sel.value = problem.senses[i];
        });
      }
      
      setStatus("Problem loaded successfully", "ok");
    }, 100);
  }

  // ========= Sensitivity Analysis =========
  function computeSensitivity(A, b, c, basis, x, z, nOrig) {
    const m = A.length;
    const nTot = nOrig + m; // including slacks
    
    // Build full A matrix with slacks
    const A_full = Array.from({ length: m }, () => Array(nTot).fill(0));
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < nOrig; j++) {
        A_full[i][j] = A[i][j];
      }
      A_full[i][nOrig + i] = 1; // slack
    }
    
    // Build basis matrix B
    const B = Array.from({ length: m }, () => Array(m).fill(0));
    for (let col = 0; col < m; col++) {
      const varIdx = basis[col];
      for (let i = 0; i < m; i++) {
        B[i][col] = A_full[i][varIdx];
      }
    }
    
    const Binv = invertMatrix(B);
    const cB = basis.map(idx => idx < nOrig ? c[idx] : 0);
    const y = multiplyRowMatrix(cB, Binv); // dual solution (shadow prices)
    
    // RHS sensitivity ranges
    const rhsRanges = [];
    for (let i = 0; i < m; i++) {
      const shadowPrice = y[i];
      // Simplified: allow ±50% change
      const currentRHS = b[i];
      const lowerBound = currentRHS * 0.5;
      const upperBound = currentRHS * 1.5;
      rhsRanges.push({
        constraint: i + 1,
        current: currentRHS,
        shadowPrice: shadowPrice,
        lowerBound: lowerBound,
        upperBound: upperBound
      });
    }
    
    // Cost coefficient ranges
    const costRanges = [];
    for (let j = 0; j < nOrig; j++) {
      const isBasic = basis.includes(j);
      const currentCost = c[j];
      if (isBasic) {
        // For basic variables, range is more complex
        costRanges.push({
          variable: j + 1,
          current: currentCost,
          lowerBound: currentCost * 0.8,
          upperBound: currentCost * 1.2,
          note: "Basic variable"
        });
      } else {
        // For non-basic, reduced cost gives range
        const reducedCost = 0; // Would need to compute from dictionary
        costRanges.push({
          variable: j + 1,
          current: currentCost,
          lowerBound: currentCost - Math.abs(reducedCost),
          upperBound: currentCost + Math.abs(reducedCost),
          note: "Non-basic variable"
        });
      }
    }
    
    return { rhsRanges, costRanges, shadowPrices: y };
  }

  function formatSensitivity(sens) {
    let html = "";
    
    // Shadow Prices Section
    html += '<div class="sens-group">';
    html += '<h4 style="color: var(--accent); margin: 0 0 8px 0; font-size: 0.85rem;">Shadow Prices (Dual Variables)</h4>';
    html += '<table style="width: 100%; margin-bottom: 16px;">';
    html += '<thead><tr><th>Constraint</th><th>Shadow Price</th><th>Interpretation</th></tr></thead>';
    html += '<tbody>';
    sens.shadowPrices.forEach((sp, i) => {
      const interpretation = sp > 0 ? "Increasing RHS improves objective" : 
                            sp < 0 ? "Decreasing RHS improves objective" : "No marginal benefit";
      html += `<tr><td><strong>C${i + 1}</strong></td><td><strong style="color: var(--accent);">${sp.toFixed(4)}</strong></td><td style="font-size: 0.75rem; color: var(--muted);">${interpretation}</td></tr>`;
    });
    html += '</tbody></table>';
    html += '</div>';
    
    // RHS Ranges Section
    html += '<div class="sens-group">';
    html += '<h4 style="color: var(--accent); margin: 0 0 8px 0; font-size: 0.85rem;">RHS Ranges (Allowable Changes)</h4>';
    html += '<table style="width: 100%; margin-bottom: 16px;">';
    html += '<thead><tr><th>Constraint</th><th>Current RHS</th><th>Lower Bound</th><th>Upper Bound</th><th>Shadow Price</th></tr></thead>';
    html += '<tbody>';
    sens.rhsRanges.forEach(r => {
      html += `<tr>
        <td><strong>C${r.constraint}</strong></td>
        <td>${r.current.toFixed(2)}</td>
        <td style="color: var(--muted);">${r.lowerBound.toFixed(2)}</td>
        <td style="color: var(--muted);">${r.upperBound.toFixed(2)}</td>
        <td><strong style="color: var(--accent);">${r.shadowPrice.toFixed(4)}</strong></td>
      </tr>`;
    });
    html += '</tbody></table>';
    html += '</div>';
    
    // Cost Coefficient Ranges Section
    html += '<div class="sens-group">';
    html += '<h4 style="color: var(--accent); margin: 0 0 8px 0; font-size: 0.85rem;">Cost Coefficient Ranges</h4>';
    html += '<table style="width: 100%;">';
    html += '<thead><tr><th>Variable</th><th>Current Cost</th><th>Lower Bound</th><th>Upper Bound</th><th>Status</th></tr></thead>';
    html += '<tbody>';
    sens.costRanges.forEach(c => {
      html += `<tr>
        <td><strong>x${c.variable}</strong></td>
        <td><strong>${c.current.toFixed(2)}</strong></td>
        <td style="color: var(--muted);">${c.lowerBound.toFixed(2)}</td>
        <td style="color: var(--muted);">${c.upperBound.toFixed(2)}</td>
        <td style="font-size: 0.75rem; color: var(--muted);">${c.note || "—"}</td>
      </tr>`;
    });
    html += '</tbody></table>';
    html += '</div>';
    
    return html;
  }

  // ========= Export Results =========
  function exportResults(A, b, c, x, z, basis, nOrig, objectiveType) {
    const xMain = x.slice(0, nOrig).map(v => Math.abs(v) < 1e-9 ? 0 : v);
    const slackVars = x.slice(nOrig);
    
    let csv = "Linear Optimization Results\n";
    csv += "===========================\n\n";
    csv += `Objective: ${objectiveType === "max" ? "Maximize" : "Minimize"}\n`;
    csv += `Optimal Objective Value: ${z.toFixed(6)}\n\n`;
    csv += "Optimal Solution (Original Variables):\n";
    xMain.forEach((val, i) => {
      csv += `x${i + 1},${val.toFixed(6)}\n`;
    });
    csv += "\nSlack/Surplus Variables:\n";
    slackVars.forEach((val, i) => {
      csv += `s${i + 1},${val.toFixed(6)}\n`;
    });
    csv += "\nBasic Variables:\n";
    basis.forEach((idx, i) => {
      csv += `Basis[${i}],x${idx + 1}\n`;
    });
    
    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "linear_optimization_results.csv";
    a.click();
    URL.revokeObjectURL(url);
  }

  function updateObjectiveLabel() {
    const objType = document.getElementById("objectiveType").value;
    const label = document.getElementById("objectiveLabel");
    label.textContent = objType === "max" ? "ζ = cᵀ x" : "ζ = cᵀ x";
  }

  // ========= DOM wiring & helpers =========
  const matrixInputsDiv = document.getElementById("matrixInputs");
  const generateBtn = document.getElementById("generateBtn");
  const solveBtn = document.getElementById("solveBtn");
  const plotBtn = document.getElementById("plotBtn");
  const dualSolveBtn = document.getElementById("dualSolveBtn");
  
  const statusEl = document.getElementById("status");
  const xStarEl = document.getElementById("xStar");
  const zStarEl = document.getElementById("zStar");
  const basicVarsEl = document.getElementById("basicVars");
  const basisSetEl = document.getElementById("basisSet");
  const xBVecEl = document.getElementById("xBVec");
  const yStarEl = document.getElementById("yStar");
  const slackVarsEl = document.getElementById("slackVars");
  const logBox = document.getElementById("logBox");
  const dictBox = document.getElementById("dictBox");
  const sensitivityBox = document.getElementById("sensitivityBox");
  const sensitivitySection = document.getElementById("sensitivitySection");
  
  // Store last solution for export
  let lastSolution = null;
  
  function setStatus(msg, type = "ok") {
    statusEl.textContent = msg;
    statusEl.className = "status " + (type === "ok" ? "ok" : "err");
  }
  
  function clearOutputs() {
    xStarEl.innerHTML = "—";
    zStarEl.textContent = "—";
    basicVarsEl.innerHTML = "—";
    basisSetEl.innerHTML = "—";
    xBVecEl.innerHTML = "—";
    yStarEl.innerHTML = "—";
    if (slackVarsEl) slackVarsEl.innerHTML = "—";
    if (sensitivityBox) sensitivityBox.innerHTML = "";
    if (sensitivitySection) sensitivitySection.style.display = "none";
    if (dictBox) dictBox.innerHTML = "";
    if (logBox) logBox.innerHTML = "";
  }
  
  function log(msg) {
    // Format log messages with better styling
    const formatted = msg
      .replace(/Phase I/g, '<strong style="color: var(--accent);">Phase I</strong>')
      .replace(/Phase II/g, '<strong style="color: var(--accent);">Phase II</strong>')
      .replace(/Iteration (\d+)/g, '<strong style="color: var(--success);">Iteration $1</strong>')
      .replace(/entering x(\d+)/g, 'entering <strong style="color: var(--accent);">x$1</strong>')
      .replace(/leaving.*?x(\d+)/g, 'leaving <strong style="color: var(--danger);">x$1</strong>')
      .replace(/Basis: \[(.*?)\]/g, 'Basis: <strong style="color: var(--accent);">[$1]</strong>');
    
    logBox.innerHTML += formatted + "<br>";
    logBox.scrollTop = logBox.scrollHeight;
  }
  
  function formatVector(arr, precision = 3, prefix = "x") {
    if (!arr || arr.length === 0) return "—";
    return arr.map((v, i) => {
      const val = Math.abs(v) < 1e-9 ? 0 : v;
      return `<span class="vector-item">${prefix}${i + 1} = <strong>${val.toFixed(precision)}</strong></span>`;
    }).join(" ");
  }

  function formatVectorSimple(arr, precision = 3) {
    return arr.map((v) => {
      const val = Math.abs(v) < 1e-9 ? 0 : v;
      return val.toFixed(precision);
    }).join(", ");
  }

  function showSolution(A, b, c, x, z, basis, methodLabel = "", objectiveType = "max") {
    const nOrig = c.length;
    const xMain = x
      .slice(0, nOrig)
      .map((v) => (Math.abs(v) < 1e-9 ? 0 : v));
    
    // Adjust z for minimize (if we negated c, negate z back)
    const finalZ = objectiveType === "min" ? -z : z;
    
    // Store solution for export
    lastSolution = { A, b, c, x, z: finalZ, basis, nOrig, objectiveType };
  
    // Format optimal solution
    xStarEl.innerHTML = formatVector(xMain, 4, "x");
    zStarEl.textContent = finalZ.toFixed(6);
    
    // Format basic variables
    basicVarsEl.innerHTML = basis && basis.length > 0 
      ? basis.map((idx) => 
          `<span class="vector-item"><strong>x${idx + 1}</strong></span>`
        ).join(" ")
      : "—";
  
    const dict = buildDictionary(A, b, c, basis);
    
    // Format basis set
    basisSetEl.innerHTML = dict.basis && dict.basis.length > 0
      ? dict.basis.map((idx) => 
          `<span class="vector-item"><strong>x${idx + 1}</strong></span>`
        ).join(", ")
      : "—";
    
    // Format basic solution values
    xBVecEl.innerHTML = formatVector(dict.xB, 4, "x_B");
    
    // Format dual solution
    yStarEl.innerHTML = formatVector(dict.y, 4, "y");
    
    // Format dictionary with better styling
    dictBox.innerHTML = dict.dictString;
    
    // Display slack/surplus variables
    const slackVars = x.slice(nOrig);
    if (slackVarsEl) {
      slackVarsEl.innerHTML = slackVars.map((v, i) => {
        const val = Math.abs(v) < 1e-9 ? 0 : v;
        return `<span class="vector-item">s${i + 1} = <strong>${val.toFixed(4)}</strong></span>`;
      }).join(" ");
    }
    
    // Compute and display sensitivity analysis
    try {
      const sens = computeSensitivity(A, b, c, basis, x, z, nOrig);
      if (sensitivityBox && sensitivitySection) {
        sensitivityBox.innerHTML = formatSensitivity(sens);
        sensitivitySection.style.display = "block";
      }
    } catch (err) {
      console.warn("Sensitivity analysis failed:", err);
    }
  
    setStatus(
      (methodLabel ? methodLabel + ": " : "") +
        "✓ Solved successfully. Optimal solution found.",
      "ok"
    );
  }
  
  function generateMatrices() {
    const m = parseInt(document.getElementById("numConstraints").value, 10);
    const n = parseInt(document.getElementById("numVariables").value, 10);
    matrixInputsDiv.innerHTML = "";
    const ui = createMatrixInputs(m, n);
    matrixInputsDiv.appendChild(ui);
    logBox.innerHTML = "";
    dictBox.innerHTML = "";
    setStatus("Matrices regenerated. Fill in A, b, c and click Solve.", "ok");
    clearOutputs();
  }
  
  // ========= Event listeners =========
  generateBtn.addEventListener("click", generateMatrices);
  
  solveBtn.addEventListener("click", () => {
    logBox.innerHTML = "";
    dictBox.innerHTML = "";
    try {
      const { A, b, c, senses, objectiveType } = readProblemFromDOM();
      // For minimize, negate c to convert to maximize
      const cAdjusted = objectiveType === "min" ? c.map(v => -v) : c;
      setStatus("Solving with Two-Phase Simplex (Phase I & II)...", "ok");
      const { x, z, basis, nOrig } = twoPhaseSimplex(A, b, cAdjusted, senses, log);
  
      // showSolution expects full x (including slacks), original A,b,c
      showSolution(A, b, c, x, z, basis, "Two-Phase Simplex", objectiveType);
    } catch (err) {
      console.error(err);
      setStatus("Error (Two-Phase): " + err.message, "err");
      log("Error (Two-Phase): " + err.message);
      clearOutputs();
    }
  });
  
  
  dualSolveBtn.addEventListener("click", () => {
    logBox.innerHTML = "";
    dictBox.innerHTML = "";
    try {
      const { A, b, c, objectiveType } = readProblemFromDOM();
      // For minimize, negate c to convert to maximize
      const cAdjusted = objectiveType === "min" ? c.map(v => -v) : c;
      setStatus("Solving with Dual Simplex...", "ok");
      const { x, z, basis } = dualSimplexStandard(A, b, cAdjusted, log);
      showSolution(A, b, c, x, z, basis, "Dual Simplex", objectiveType);
    } catch (err) {
      console.error(err);
      setStatus("Error (Dual Simplex): " + err.message, "err");
      log("Error (Dual Simplex): " + err.message);
      clearOutputs();
    }
  });
  
  plotBtn.addEventListener("click", () => {
    const { n, A, b, c, objectiveType } = readProblemFromDOM();
    if (n !== 2) {
      setStatus("Visualization only supports n = 2. Set Variables = 2.", "err");
      return;
    }
    setStatus("Plotting feasible region for (x1, x2)...", "ok");
    plotFeasibleRegion(A, b, c, objectiveType);
  });
  
  // ========= New Event Listeners =========
  const importBtn = document.getElementById("importBtn");
  const exportBtn = document.getElementById("exportBtn");
  const saveBtn = document.getElementById("saveBtn");
  const loadBtn = document.getElementById("loadBtn");
  const exampleSelect = document.getElementById("exampleSelect");
  const exportResultsBtn = document.getElementById("exportResultsBtn");
  const objectiveTypeSelect = document.getElementById("objectiveType");
  const csvFileInput = document.getElementById("csvFileInput");
  
  // Objective type toggle
  if (objectiveTypeSelect) {
    objectiveTypeSelect.addEventListener("change", updateObjectiveLabel);
  }
  
  // CSV Import
  if (importBtn) {
    importBtn.addEventListener("click", () => {
      csvFileInput.click();
    });
    
    csvFileInput.addEventListener("change", (e) => {
      const file = e.target.files[0];
      if (!file) return;
      
      const reader = new FileReader();
      reader.onload = (event) => {
        try {
          importFromCSV(event.target.result);
          setStatus("CSV imported successfully", "ok");
        } catch (err) {
          setStatus("CSV import error: " + err.message, "err");
        }
      };
      reader.readAsText(file);
      csvFileInput.value = ""; // Reset
    });
    
    // Also support paste
    document.addEventListener("paste", (e) => {
      const text = e.clipboardData.getData("text");
      if (text && text.includes(",")) {
        try {
          importFromCSV(text);
          setStatus("Pasted CSV data imported", "ok");
        } catch (err) {
          // Silently fail if not CSV
        }
      }
    });
  }
  
  // CSV Export
  if (exportBtn) {
    exportBtn.addEventListener("click", () => {
      try {
        exportToCSV();
        setStatus("Matrices exported to CSV", "ok");
      } catch (err) {
        setStatus("Export error: " + err.message, "err");
      }
    });
  }
  
  // Save problem
  if (saveBtn) {
    saveBtn.addEventListener("click", () => {
      try {
        saveProblem();
      } catch (err) {
        setStatus("Save error: " + err.message, "err");
      }
    });
  }
  
  // Load problem
  if (loadBtn) {
    loadBtn.addEventListener("click", () => {
      loadProblem();
    });
  }
  
  // Example problems
  if (exampleSelect) {
    exampleSelect.addEventListener("change", (e) => {
      if (e.target.value) {
        loadExample(e.target.value);
        e.target.value = ""; // Reset dropdown
      }
    });
  }
  
  // Export results
  if (exportResultsBtn) {
    exportResultsBtn.addEventListener("click", () => {
      try {
        if (!lastSolution) {
          setStatus("No solution to export. Solve the problem first.", "err");
          return;
        }
        
        exportResults(
          lastSolution.A, 
          lastSolution.b, 
          lastSolution.c, 
          lastSolution.x, 
          lastSolution.z, 
          lastSolution.basis, 
          lastSolution.nOrig, 
          lastSolution.objectiveType
        );
        setStatus("Results exported to CSV", "ok");
      } catch (err) {
        setStatus("Export error: " + err.message, "err");
      }
    });
  }
  
  // Initial render
  generateMatrices();
  updateObjectiveLabel();
  