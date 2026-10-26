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
function widget() {
  const handlers = {};
  return {handlers, on: (key, fn) => { handlers[key] = fn; }};
}
const flush = () => new Promise(resolve => setImmediate(resolve));
function selectModel(el, name) {
  const frame = el._transitionData._frames.find(frame => frame.name === name);
  el.data = frame.data;
  Object.assign(el.layout, frame.layout);
  el.handlers.plotly_animatingframe({name});
}
(async () => {
  const a = widget(), b = widget();
  draw(a, {layout:{}}, JSON.parse(JSON.stringify(source)));
  draw(b, {layout:{}}, JSON.parse(JSON.stringify(source)));
  await flush();
  assert(a.data.some(trace => trace.x.length));
  for (let component = 0; component < 3; component++) {
    a.handlers.plotly_buttonclicked({button:{method:'skip',args:[component]}});
    await flush();
    assert.equal(a.layout.sliders[0].steps.length, 3);
    assert.equal(a.layout.updatemenus.length, 1);
    assert(!JSON.stringify(a.layout.updatemenus).includes('Play'));
    assert(!JSON.stringify(a.layout.updatemenus).includes('Pause'));
    const beforeAnimation = banks;
    for (const model of source.models) {
      selectModel(a, model);
      const available = model !== 'scaa';
      assert.equal(a.data.some(trace => trace.x.length), available);
      assert.equal(a.layout.annotations.length, available ? 0 : 1);
      if (!available) assert.equal(a.layout.annotations[0].text, 'Unavailable for this model component');
      a.data.forEach(trace => assert(Math.abs(trace.marker.sizeref - 4/900) < 1e-15));
      if (available && component === 2) assert.deepEqual(a.layout.yaxis.ticktext,
        model === 'sam' ? ['2-4','5-8'] : ['2','3-5','6-8']);
    }
    assert.equal(banks, beforeAnimation); // Slider selection must not rebuild the frame bank.
    a.handlers.plotly_buttonclicked({button:{method:'skip',args:[(component+1)%3]}});
    await flush();
    assert.equal(a.layout.sliders[0].active, 2); // retain the selected model
    assert.equal(b.layout.sliders[0].active, 1); // independent widgets
  }
  for (let component = 0; component < 3; component++)
    a.handlers.plotly_buttonclicked({button:{method:'skip',args:[component]}});
  await flush(); await flush();
  assert.equal(a.layout.updatemenus[0].active, 2);
  // htmlwidgets may unbox one-element coordinates and state-block labels.
  const singleConfig = JSON.parse(JSON.stringify(source));
  singleConfig.models = 'single'; singleConfig.model = 0; singleConfig.component = 2;
  singleConfig.components.forEach(component => {
    component.cells = [component.cells[1]];
    component.cells[0].layout.yaxis.tickvals = 1;
    component.cells[0].layout.yaxis.ticktext = '2-14';
    component.cells[0].data.forEach(trace => {
      trace.x = 2000; trace.y = 1; trace.text = 'one'; trace.marker.size = 1;
    });
  });
  const single = widget(); draw(single, {layout:{}}, singleConfig); await flush();
  assert.deepEqual(single.layout.yaxis.ticktext, ['2-14']);
  assert.deepEqual(single.data[0].x, [2000]);
  assert.deepEqual(single.layout.sliders, []);
  console.log('Process bubble controls passed');
})().catch(error => { console.error(error); process.exit(1); });
