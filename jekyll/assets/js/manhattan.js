function processData(rows) {
  const toNumber = value => {
    const number = Number(value);
    return value && Number.isFinite(number) ? number : null;
  };
  const positions = rows.map(row => toNumber(row.POS));
  const snps = rows.map(row => row.SNP);
  const marginal = rows.map(row => toNumber(row["GWAS.LOGP"]));
  const conditional = rows.map(row => toNumber(row["GWAS_cond.LOGP"]));
  const maxLogP = marginal.concat(conditional).reduce(
    (max, value) => value === null ? max : Math.max(max, value), 0
  );
  const range = [0, Math.max(1, maxLogP * 1.05)];
  const markers = {
    x: positions,
    text: snps,
    mode: "markers",
    hovertemplate: "%{text}<br>Position %{x}<br>−log₁₀ P %{y}<extra></extra>",
    marker: { size: 5, color: "#333" },
  };

  return Plotly.newPlot("graph", [
    { ...markers, y: marginal, name: "Marginal" },
    { ...markers, y: conditional, name: "Conditional", xaxis: "x2", yaxis: "y2" },
  ], {
    height: 520,
    margin: { l: 65, r: 25, b: 55, t: 30 },
    xaxis: { anchor: "y", matches: "x2", showticklabels: false },
    xaxis2: { anchor: "y2", title: { text: "Physical position" } },
    yaxis: { domain: [0.58, 1], range, matches: "y2", title: { text: "−log₁₀ P" } },
    yaxis2: { domain: [0, 0.42], range, title: { text: "−log₁₀ P" } },
    annotations: [
      { xref: "paper", yref: "paper", x: 0, y: 1, yanchor: "bottom", showarrow: false, text: "Marginal GWAS" },
      { xref: "paper", yref: "paper", x: 0, y: 0.42, yanchor: "bottom", showarrow: false, text: "Conditional GWAS" },
    ],
    hovermode: "closest",
    showlegend: false,
  }, { responsive: true });
}
