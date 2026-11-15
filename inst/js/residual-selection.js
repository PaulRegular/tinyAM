function(el, x, config) {
  'use strict';
  function array(value) { return Array.isArray(value) ? value : [value]; }
  var frames = array(config.frames);
  var defaults = array(config.defaultViews);
  var currentFrame = 0;
  var currentView = config.view;
  var automatic = config.automatic;
  var drawing = false;
  var pending = false;
  var controls = document.createElement('div');
  controls.className = 'tam-residual-controls';
  controls.style.cssText = 'position:absolute;top:26px;left:8px;right:8px;z-index:5;' +
    'display:flex;gap:.8em;align-items:center;pointer-events:none;';

  function select(label, name, values, labels, selected, change) {
    var field = document.createElement('label');
    field.style.cssText = 'display:flex;gap:.4em;align-items:center;min-width:0;' +
      'white-space:nowrap;pointer-events:auto;margin:0;';
    var description = document.createElement('span');
    description.textContent = label;
    var input = document.createElement('select');
    input.className = 'tam-residual-' + name;
    input.style.cssText = 'min-width:0;font:inherit;';
    values.forEach(function(value, i) {
      var option = document.createElement('option');
      option.value = String(value);
      option.textContent = labels[i];
      input.appendChild(option);
    });
    input.value = String(selected);
    input.addEventListener('change', function() {
      change(Number(input.value));
      draw();
    });
    field.appendChild(description);
    field.appendChild(input);
    controls.appendChild(field);
    return input;
  }

  var typeSelect;
  if (config.showSelector) {
    var types = array(config.selectors);
    typeSelect = select('Residuals', 'type', types,
      types.map(function(i) { return config.labels[i]; }), currentView, function(view) {
        currentView = view;
        automatic = false;
      });
  }
  if (frames.length > 1) {
    select(config.frameLabel, 'frame', frames.map(function(_, i) { return i; }),
      frames, currentFrame, function(frame) {
        currentFrame = frame;
        if (automatic) {
          currentView = defaults[frame];
          if (typeSelect) typeSelect.value = String(currentView);
        }
      });
  }
  if (controls.childElementCount) {
    if (window.getComputedStyle(el).position === 'static') el.style.position = 'relative';
    el.appendChild(controls);
  }

  function draw() {
    pending = true;
    if (drawing) return;
    drawing = true;
    pending = false;
    // Plotly can update input objects. Fresh copies keep other views independent.
    var cell = JSON.parse(JSON.stringify(config.views[currentView].cells[currentFrame]));
    var data = cell.data || [];
    data.forEach(function(trace) {
      // htmlwidgets serializes length-one vectors as scalars.
      ['x', 'y', 'text', 'hovertext', 'customdata', 'ids', 'key'].forEach(function(key) {
        if (trace[key] !== undefined) trace[key] = array(trace[key]);
      });
      if (trace.marker && trace.marker.size !== undefined) trace.marker.size = array(trace.marker.size);
    });
    var layout = Object.assign({}, x.layout, cell.layout, {
      autosize: true, margin: Object.assign({}, cell.layout.margin,
        {t: controls.childElementCount ? 85 : 60})
    });
    function finished() {
      drawing = false;
      if (pending) draw();
    }
    Promise.resolve().then(function() { return Plotly.react(el, data, layout); })
      .then(finished, function(error) { console.warn(error); finished(); });
  }
  draw();
}
