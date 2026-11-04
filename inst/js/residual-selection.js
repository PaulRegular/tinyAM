function(el, x, config) {
  'use strict';
  function array(value) { return Array.isArray(value) ? value : [value]; }
  config.frames = array(config.frames);
  config.defaultViews = array(config.defaultViews);
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
  var state = {view: config.view, frame: config.frame, automatic: config.automatic}, busy = false, pending = false;
  function selectedView(frame) {
    return state.automatic ? config.defaultViews[frame] : state.view;
  }
  function draw() {
    if (busy) { pending = true; return; }
    busy = true;
    state.view = selectedView(state.frame);
    var cell = config.views[state.view].cells[state.frame];
    var frames = config.frames.map(function(name, i) {
      var frameCell = config.views[selectedView(i)].cells[i];
      return {name: name, data: frameCell.data, layout: frameCell.layout,
        traces: frameCell.data.map(function(trace, j) { return j; })};
    });
    var layout = Object.assign({}, x.layout, cell.layout, {
      autosize: true, margin: Object.assign({}, cell.layout.margin, {t: config.showSelector ? 100 : 60}),
      // Leave the upper-right modebar clear and reserve a separate method row.
      updatemenus: config.showSelector ? [{name: 'tam-residual-type', type: 'dropdown',
        x: 0, y: 1, xanchor: 'left', yanchor: 'bottom', pad: {b: 28}, active: state.view,
        buttons: ['Standardized', 'OSA'].map(function(label, i) {
          return {label: label, method: 'skip', execute: false, args: [i]};
        })}] : [],
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
    state.automatic = false;
    draw();
  });
  el.on('plotly_animatingframe', function(event) {
    var index = config.frames.indexOf(event.name);
    if (index >= 0) {
      state.frame = index;
      if (state.automatic) {
        state.view = config.defaultViews[index];
        if (!busy && config.showSelector) Plotly.relayout(el, {'updatemenus[0].active': state.view});
      }
      if (busy) pending = true;
    }
  });
  draw();
}
