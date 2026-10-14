function(el, x, config) {
  // htmlwidgets unboxes length-one vectors when serializing to JSON.
  function array(value) { return Array.isArray(value) ? value : [value]; }
  config.models = array(config.models);
  config.metrics.forEach(function(metric) {
    metric.years = array(metric.years);
    metric.ages = array(metric.ages);
    metric.cells.forEach(function(cell) {
      ['linear', 'log', 'text'].forEach(function(field) {
        cell[field] = array(cell[field]).map(array);
      });
    });
  });
  var state = {metric: 0, model: 0, scale: 'linear', type: config.type};
  var camera = {eye: {x: 2, y: 2, z: 1.6}}, busy = false, pending = false;
  function menu(name, labels, values, active, top) {
    return {name: name, type: 'buttons', x: 0, y: 1, xanchor: 'right',
      yanchor: 'top', pad: {r: 60, t: top}, active: active,
      buttons: labels.map(function(label, i) {
        return {label: label, method: 'skip', execute: false,
          args: [name, values[i]]};
      })};
  }
  function trace(model) {
    var metric = config.metrics[state.metric], cell = metric.cells[model];
    var limits = metric[state.scale];
    var values = cell[state.scale];
    var available = values.some(function(row) {
      return row.some(function(value) { return value !== null && isFinite(value); });
    });
    return {type: state.type, x: metric.years, y: metric.ages,
      z: values, visible: available, text: cell.text, hovertemplate: '%{text}<extra></extra>',
      colorscale: 'Viridis', zmin: limits[0], zmax: limits[1],
      cmin: limits[0], cmax: limits[1], connectgaps: false,
      colorbar: {title: {text: (state.scale === 'log' ? 'log10 ' : '') +
        metric.name + (metric.unit ? ' (' + metric.unit + ')' : '')}},
      contours: state.type === 'contour' ? {coloring: 'heatmap'} : undefined,
      showscale: true};
  }
  function draw() {
    if (busy) { pending = true; return; }
    busy = true;
    var metric = config.metrics[state.metric], limits = metric[state.scale];
    var yearRange = [metric.years[0] - .5, metric.years[metric.years.length - 1] + .5];
    var ageRange = [metric.ages[0] - .5, metric.ages[metric.ages.length - 1] + .5];
    var menus = [menu('scale', ['Linear', 'Log'], ['linear', 'log'],
      state.scale === 'linear' ? 0 : 1, 0),
      menu('type', ['Heatmap', 'Surface', 'Contour'], ['heatmap', 'surface', 'contour'],
        ['heatmap', 'surface', 'contour'].indexOf(state.type), 75)];
    if (config.metrics.length > 1) menus.push(menu('metric',
      config.metrics.map(function(m) { return m.name; }),
      config.metrics.map(function(m, i) { return i; }), state.metric, 180));
    var frames = config.models.map(function(model, i) {
      return {name: model, data: [trace(i)], traces: [0]};
    });
    var data = [trace(state.model)];
    var layout = Object.assign({}, x.layout, {autosize: true, margin: {l: 170, r: 85, t: 25, b: 100},
      uirevision: 'matrix', updatemenus: menus, xaxis: {title: 'Year', range: yearRange},
      yaxis: {title: 'Age', range: ageRange},
      scene: {xaxis: {title: 'Year', range: yearRange}, yaxis: {title: 'Age', range: ageRange},
        zaxis: {title: (state.scale === 'log' ? 'log10 ' : '') + metric.name,
          range: limits}, camera: camera, uirevision: 'matrix-camera'},
      sliders: config.models.length > 1 ? [{active: state.model,
        currentvalue: {prefix: 'Model: '}, pad: {t: 45},
        steps: config.models.map(function(model) {
          return {label: model, method: 'animate', args: [[model],
            {mode: 'immediate', frame: {duration: 0, redraw: true},
              transition: {duration: 0}}]};
        })}] : []});
    var oldFrames = el._transitionData && el._transitionData._frames;
    var remove = oldFrames && oldFrames.length ?
      Plotly.deleteFrames(el) : Promise.resolve();
    remove.then(function() { return Plotly.react(el, data, layout); })
      .then(function() { return Plotly.addFrames(el, frames); })
      .then(function() {
        busy = false;
        if (pending) { pending = false; draw(); }
      });
  }
  el.on('plotly_buttonclicked', function(event) {
    if (!event.button || !event.button.args) return;
    var args = event.button.args;
    if (['metric', 'scale', 'type'].indexOf(args[0]) < 0) return;
    state[args[0]] = args[1];
    draw();
  });
  el.on('plotly_animatingframe', function(event) {
    var index = config.models.indexOf(event.name);
    if (index >= 0) state.model = index;
  });
  el.on('plotly_relayout', function(event) {
    if (event['scene.camera']) camera = event['scene.camera'];
  });
  // Hidden dashboard tabs need a resize when they become visible.
  var visible = false;
  var observer = new IntersectionObserver(function(entries) {
    entries.forEach(function(entry) {
      if (entry.isIntersecting && !visible) Plotly.Plots.resize(el);
      visible = entry.isIntersecting;
    });
  });
  observer.observe(el);
  draw();
}
