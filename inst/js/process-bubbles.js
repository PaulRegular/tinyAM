function(el, x, config) {
  function array(value) { return Array.isArray(value) ? value : [value]; }
  config.models = array(config.models);
  config.components.forEach(function(component) {
    component.cells.forEach(function(cell) {
      ['tickvals', 'ticktext'].forEach(function(key) {
        cell.layout.yaxis[key] = array(cell.layout.yaxis[key]);
      });
      cell.data.forEach(function(trace) {
        ['x', 'y', 'text'].forEach(function(key) {
          if (trace[key] !== undefined) trace[key] = array(trace[key]);
        });
        trace.marker.size = array(trace.marker.size);
      });
    });
  });
  var state = {component: config.component, model: config.model};
  var busy = false, pending = false;
  function draw() {
    if (busy) { pending = true; return; }
    busy = true;
    var component = config.components[state.component];
    var cell = component.cells[state.model];
    var frames = config.models.map(function(model, i) {
      return {name: model, data: component.cells[i].data,
        layout: component.cells[i].layout, traces: [0, 1]};
    });
    var layout = Object.assign({}, x.layout, cell.layout, {
      title: '', showlegend: false, autosize: true,
      margin: {l: 150, r: 25, t: 25, b: 110},
      updatemenus: [{type: 'buttons', x: 0, y: 1, xanchor: 'right',
        yanchor: 'top', pad: {r: 45}, active: state.component,
        buttons: config.components.map(function(c, i) {
          return {label: c.name, method: 'skip', execute: false, args: [i]};
        })}],
      sliders: config.models.length > 1 ? [{active: state.model, pad: {t: 50},
        currentvalue: {prefix: 'Model: '}, steps: config.models.map(function(model) {
          return {label: model, method: 'animate', args: [[model],
            {mode: 'immediate', frame: {duration: 0, redraw: true}, transition: {duration: 0}}]};
        })}] : []
    });
    // Only the component buttons replace the frame bank. Model selection stays native.
    var pause = el._transitionData ? Plotly.animate(el, [null],
      {mode: 'immediate', frame: {duration: 0}, transition: {duration: 0}}) : Promise.resolve();
    pause.then(function() {
      var frames = el._transitionData && el._transitionData._frames;
      return frames && frames.length ? Plotly.deleteFrames(el) : null;
    }).then(function() { return Plotly.react(el, cell.data, layout); })
      .then(function() { return Plotly.addFrames(el, frames); })
      .then(function() {
        busy = false;
        if (pending) { pending = false; draw(); }
      });
  }
  el.on('plotly_buttonclicked', function(event) {
    if (!event.button || event.button.method !== 'skip') return;
    state.component = event.button.args[0];
    draw();
  });
  el.on('plotly_animatingframe', function(event) {
    var index = config.models.indexOf(event.name);
    if (index >= 0) state.model = index;
  });
  draw();
}
