(function () {
  "use strict";
  const root = document.currentScript.parentElement;
  const buttons = root.querySelectorAll(".tam-age-year-buttons button");
  const panels = root.querySelectorAll(".tam-age-year-panel");
  function resize() {
    panels.forEach(function (panel) {
      if (panel.style.display === "none") return;
      const graph = panel.querySelector(".js-plotly-plot");
      if (graph && graph._fullLayout && window.Plotly && panel.clientWidth) {
        Plotly.Plots.resize(graph);
      }
    });
  }
  buttons.forEach(function (button) {
    button.addEventListener("click", function () {
      panels.forEach(function (panel) {
        panel.style.display = panel.dataset.view === button.dataset.view ? "" : "none";
      });
      buttons.forEach(function (other) {
        const active = other === button;
        other.setAttribute("aria-pressed", String(active));
        other.classList.toggle("active", active);
      });
      requestAnimationFrame(resize);
    });
    button.classList.toggle("active", button.getAttribute("aria-pressed") === "true");
  });
  // Dashboard tabs can reveal a plot without a window resize.
  const observer = new ResizeObserver(function () { requestAnimationFrame(resize); });
  observer.observe(root);
})();
