function(el, x, config) {
  'use strict';
  function array(value) { return Array.isArray(value) ? value : [value]; }
  config.frames = array(config.frames);
  config.views.forEach(function(view) {
    view.cells.forEach(function(cell) {
      cell.data.forEach(function(trace) {
        ['x', 'y', 'text', 'hovertext', 'customdata'].forEach(function(key) {
          if (trace[key] !== undefined) trace[key] = array(trace[key]);
        });
        if (trace.marker && trace.marker.size !== undefined) trace.marker.size = array(trace.marker.size);
      });
    });
  });
  var state = {view: config.view, frame: config.frame}, busy = false, pending = false;
  function draw() {
    if (busy) { pending = true; return; }
    busy = true;
    var view = config.views[state.view], cell = view.cells[state.frame];
    var frames = config.frames.map(function(name, i) {
      return {name: name, data: view.cells[i].data, layout: view.cells[i].layout,
        traces: view.cells[i].data.map(function(trace, j) { return j; })};
    });
    var layout = Object.assign({}, x.layout, cell.layout, {
      autosize: true, margin: Object.assign({}, cell.layout.margin, {t: 100}),
      // Leave the upper-right modebar clear and reserve a separate method row.
      updatemenus: [{name: 'tam-residual-type', type: 'dropdown',
        x: 0, y: 1, xanchor: 'left', yanchor: 'bottom', pad: {b: 28}, active: state.view,
        buttons: ['Auto', 'Standardized', 'OSA'].map(function(label, i) {
          return {label: label, method: 'skip', execute: false, args: [i]};
        })}],
      sliders: config.animated ? [{active: state.frame, pad: {t: 45},
        currentvalue: {prefix: config.prefix}, steps: config.frames.map(function(name) {
          return {label: name, method: 'animate', args: [[name],
            {mode: 'immediate', frame: {duration: 0, redraw: true}, transition: {duration: 0}}]};
        })}] : []
    });
    var pause = el._transitionData ? Plotly.animate(el, [null],
      {mode: 'immediate', frame: {duration: 0}, transition: {duration: 0}}) : Promise.resolve();
    pause.then(function() {
      var previous = el._transitionData && el._transitionData._frames;
      return previous && previous.length ? Plotly.deleteFrames(el) : null;
    }).then(function() { return Plotly.react(el, cell.data, layout); })
      .then(function() { return config.animated ? Plotly.addFrames(el, frames) : null; })
      .then(function() {
        busy = false;
        if (pending) { pending = false; draw(); }
      }).catch(function(error) { busy = false; console.warn(error); });
  }
  el.on('plotly_buttonclicked', function(event) {
    if (!event.menu || event.menu.name !== 'tam-residual-type') return;
    state.view = event.button.args[0];
    draw();
  });
  el.on('plotly_animatingframe', function(event) {
    var index = config.frames.indexOf(event.name);
    if (index >= 0) {
      state.frame = index;
      if (busy) pending = true;
    }
  });
  draw();
}
