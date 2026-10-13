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
        Plotly.Plots.resize(graph).catch(function (error) {
          // A view may be hidden while Plotly's deferred resize is pending.
          if (panel.clientWidth) console.warn(error);
        });
      }
    });
  }
  buttons.forEach(function (button) {
    button.addEventListener("click", function () {
      panels.forEach(function (panel) {
        if (panel.style.display !== "none" && panel.dataset.view !== button.dataset.view) {
          const graph = panel.querySelector(".js-plotly-plot");
          if (graph && graph._fullLayout && graph.querySelector(".slider-container")) {
            Plotly.animate(graph, [null], {mode: "immediate", frame: {duration: 0},
              transition: {duration: 0}}).catch(function () {
                // Cancelling an animation can reject its pending promise.
              });
          }
        }
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
