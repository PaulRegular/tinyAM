(function () {
  "use strict";
  const root = document.currentScript.parentElement;
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
  function switchView(view) {
    panels.forEach(function (panel) {
      if (panel.style.display !== "none" && panel.dataset.view !== view) {
        const graph = panel.querySelector(".js-plotly-plot");
        if (graph && graph._fullLayout && graph.querySelector(".slider-container")) {
          Plotly.animate(graph, [null], {mode: "immediate", frame: {duration: 0},
            transition: {duration: 0}}).catch(function () {
              // Cancelling an animation can reject its pending promise.
            });
        }
      }
      panel.style.display = panel.dataset.view === view ? "" : "none";
    });
    requestAnimationFrame(resize);
  }
  HTMLWidgets.addPostRenderHandler(function () {
    panels.forEach(function (panel) {
      const graph = panel.querySelector(".js-plotly-plot");
      if (graph) graph.on("plotly_buttonclicked", function (event) {
        if (event.menu.name === "tam-age-year-view") switchView(event.button.args[0]);
      });
    });
    resize();
  });
  // Dashboard tabs can reveal a plot without a window resize.
  const observer = new ResizeObserver(function () { requestAnimationFrame(resize); });
  observer.observe(root);
})();
