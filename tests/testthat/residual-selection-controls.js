const fs = require('fs'), assert = require('assert');
const draw = eval('(' + fs.readFileSync(process.argv[2], 'utf8') + ')');
const source = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
let banks = 0;
global.Plotly = {
  react: async (el, data, layout) => { el.data = data; el.layout = layout; },
  addFrames: async (el, frames) => { banks++; el._transitionData = {_frames: frames}; },
  deleteFrames: async el => { el._transitionData = {_frames: []}; },
  animate: async () => {}
};
const flush = () => new Promise(resolve => setImmediate(resolve));
function widget() {
  const handlers = {};
  return {handlers, on: (name, fun) => { handlers[name] = fun; }};
}
function chooseType(el, index) {
  el.handlers.plotly_buttonclicked({menu: {name: 'tam-residual-type'}, button: {args: [index]}});
}
function chooseFrame(el, name) {
  const frame = el._transitionData._frames.find(frame => frame.name === name);
  el.data = frame.data;
  Object.assign(el.layout, frame.layout);
  el.handlers.plotly_animatingframe({name});
}
const points = el => el.data.reduce((n, trace) => n + trace.x.length, 0);
(async () => {
  const a = widget(), b = widget();
  draw(a, {layout: {}}, JSON.parse(JSON.stringify(source)));
  draw(b, {layout: {}}, JSON.parse(JSON.stringify(source)));
  await flush();
  assert.equal(points(a), 2);
  assert.equal(points(b), 2);
  assert.deepEqual(a.layout.updatemenus[0].buttons.map(button => button.label), ['Auto', 'Standardized', 'OSA']);
  assert.equal(a.layout.updatemenus[0].x, 0);
  assert.equal(a.layout.updatemenus[0].xanchor, 'left');
  assert.equal(a.layout.updatemenus[0].yanchor, 'bottom');
  assert.equal(a.layout.updatemenus[0].pad.b, 28);
  assert.equal(a.layout.margin.t, 100); // Header rows stay clear of the upper-right modebar.
  chooseFrame(a, 'B');
  assert.equal(points(a), 6); // Auto keeps a model without cached OSA.
  const before = banks;
  chooseFrame(a, 'A');
  chooseFrame(a, 'B');
  assert.equal(banks, before); // Native frame changes keep their frame bank.
  chooseType(a, 2); await flush();
  assert.equal(a.layout.sliders[0].active, 1);
  assert.equal(points(a), 0);
  assert(a.layout.annotations[0].text.includes('OSA residuals unavailable'));
  assert.equal(points(b), 2); // Widgets have independent selections.
  chooseType(a, 1); await flush();
  assert.equal(points(a), 6);
  assert.equal(a.layout.sliders[0].active, 1);
  chooseFrame(a, 'A');
  assert.equal(points(a), 6);
  chooseType(a, 0); await flush();
  assert.equal(points(a), 2); // Partial OSA never fills its gaps from standard residuals.
  assert(a.layout.annotations[0].text.includes('oneStepGeneric'));
  for (let i = 0; i < 3; i++) chooseType(a, i);
  await flush(); await flush();
  assert.equal(a.layout.updatemenus[0].active, 2);
  const singleConfig = JSON.parse(JSON.stringify(source));
  singleConfig.frames = 'Only'; singleConfig.animated = false;
  singleConfig.views.forEach(view => {
    view.cells = [view.cells[0]];
    view.cells[0].data.forEach(trace => {
      trace.x = 2000; trace.y = 1; trace.text = 'one'; trace.customdata = 'custom';
      trace.marker = trace.marker || {}; trace.marker.size = 1;
    });
  });
  const single = widget(); draw(single, {layout: {}}, singleConfig); await flush();
  assert.deepEqual(single.data[0].x, [2000]);
  assert.deepEqual(single.data[0].customdata, ['custom']);
  assert.deepEqual(single.layout.sliders, []);
  console.log('Residual controls passed');
})().catch(error => { console.error(error); process.exit(1); });
