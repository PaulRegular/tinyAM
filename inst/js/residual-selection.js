function(el, x, config) {
  'use strict';
  function array(value) { return Array.isArray(value) ? value : [value]; }
  var frames = array(config.frames);
  var defaults = array(config.defaultViews);
  var currentFrame = config.frame;
  var currentView = config.view;
  var automatic = config.automatic;
  var drawing = false;
  var pending = false;
  var controls = document.createElement('div');
  controls.className = 'tam-residual-controls';
  controls.style.cssText = 'position:absolute;left:8px;right:8px;z-index:5;' +
    'display:flex;gap:.8em;align-items:center;pointer-events:none;';
  controls.style.top = config.compactHeader ? '8px' : '26px';

  function setOptions(input, values, labels, selected) {
    input.replaceChildren();
    values.forEach(function(value, i) {
      var option = document.createElement('option');
      option.value = String(value);
      option.textContent = labels[i];
      input.appendChild(option);
    });
    input.value = String(selected);
  }

  function select(label, name, change) {
    var field = document.createElement('label');
    field.style.cssText = 'display:flex;gap:.4em;align-items:center;min-width:0;' +
      'white-space:nowrap;pointer-events:auto;margin:0;';
    var description = document.createElement('span');
    description.textContent = label;
    var input = document.createElement('select');
    input.className = 'tam-residual-' + name;
    input.style.cssText = 'min-width:0;font:inherit;';
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
  function updateViewOptions() {
    if (!typeSelect) return;
    var options = config.viewOptions[currentFrame];
    // Multiple views can share a label while preserving the user's choice.
    setOptions(typeSelect, array(options.values), array(options.labels), array(options.display)[currentView]);
  }
  if (config.showSelector) {
    typeSelect = select(config.viewLabel, 'type', function(view) {
      currentView = view;
      automatic = false;
    });
    updateViewOptions();
  }
  if (frames.length > 1) {
    var modelSelect = select(config.frameLabel, 'frame', function(frame) {
      currentFrame = frame;
      if (automatic) currentView = defaults[frame];
      updateViewOptions();
    });
    setOptions(modelSelect, frames.map(function(_, i) { return i; }), frames, currentFrame);
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
    ['xaxis', 'yaxis'].forEach(function(axis) {
      ['tickvals', 'ticktext'].forEach(function(key) {
        if (cell.layout[axis] && cell.layout[axis][key] != null) {
          cell.layout[axis][key] = array(cell.layout[axis][key]);
        }
      });
    });
    var top = controls.childElementCount ? 85 : 60;
    if (config.compactHeader) top -= 25;
    var layout = Object.assign({}, x.layout, cell.layout, {
      autosize: true, margin: Object.assign({}, cell.layout.margin,
        {t: top})
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
