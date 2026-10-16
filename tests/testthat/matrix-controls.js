// Exercise the widget coordinator without a browser or a graphics driver.
const fs = require('fs');
const assert = require('assert');
const draw = eval('(' + fs.readFileSync(process.argv[2], 'utf8') + ')');
global.IntersectionObserver = class { observe() {} };
global.Plotly = {
  react: async (el, data, layout) => { el.data = data; el.layout = layout; },
  addFrames: async (el, frames) => { el._transitionData = {_frames: frames}; },
  deleteFrames: async el => { el._transitionData = {_frames: []}; },
  Plots: {resize: () => {}}
};
function widget() {
  const handlers = {};
  return {on: (event, fn) => { handlers[event] = fn; }, handlers};
}
const config = {type: 'heatmap', models: ['tinyAM', 'Official'],
  colorscale: [[0, '#003e56'], [1, '#f05c28']],
  metrics: ['Abundance', 'Biomass', 'SSB', 'F', 'M', 'Z'].map((name, i) => ({
    name, unit: i > 2 ? 'per year' : i ? 'kg' : 'fish', years: [2000, 2001], ages: [1, 2],
    linear: [0, 100], log: [0, 2], cells: [0, 1].map(model => ({
      linear: [[0, 1 + model], [null, 100]], log: [[null, model], [null, 2]],
      text: [['zero', 'natural value'], ['missing', '100 kg']]
    }))
  }))};
const flush = () => new Promise(resolve => setImmediate(resolve));
(async () => {
  config.metrics[2].cells[1].linear = [[null, null], [null, null]];
  config.metrics[2].cells[1].log = [[null, null], [null, null]];
  config.metrics[5].cells[1].linear = [[null, null], [null, null]];
  config.metrics[5].cells[1].log = [[null, null], [null, null]];
  const a = widget(), b = widget();
  draw(a, {}, config); draw(b, {}, config); await flush();
  a.handlers.plotly_animatingframe({name: 'Official'});
  const camera = {eye: {x: 1, y: 2, z: 3}};
  a.handlers.plotly_relayout({'scene.camera': camera});
  for (let metric = 0; metric < config.metrics.length; metric++) {
    for (const scale of ['linear', 'log']) for (const type of ['heatmap', 'surface', 'contour']) {
      for (const [key, value] of Object.entries({metric, scale, type})) {
        a.handlers.plotly_buttonclicked({button: {args: [key, value]}});
        await flush();
      }
      assert.equal(a.data[0].type, type);
      assert.deepEqual(a.data[0].colorscale, config.colorscale);
      assert.deepEqual(a.data[0].z, config.metrics[metric].cells[1][scale]);
      assert.equal(a.data[0].visible, metric !== 2 && metric !== 5);
      assert.deepEqual(a.layout.xaxis.range, [1999.5, 2001.5]);
      assert.deepEqual(a.layout.yaxis.range, [.5, 2.5]);
      assert.equal(a.layout.sliders[0].active, 1);
      assert.deepEqual(a.layout.scene.camera, camera);
      assert.deepEqual(a.layout.scene.zaxis.range, config.metrics[metric][scale]);
      assert.equal(a._transitionData._frames.length, 2);
      assert.deepEqual(a._transitionData._frames[0].data[0].z,
        config.metrics[metric].cells[0][scale]);
      assert.deepEqual(a._transitionData._frames[0].data[0].colorscale, config.colorscale);
      assert.equal(b.data[0].type, 'heatmap');
      assert.equal(b.layout.sliders[0].active, 0);
    }
  }
  // Rapid clicks during a redraw must retain the last requested combination.
  for (const [key, value] of Object.entries({metric: 0, scale: 'linear', type: 'surface'})) {
    a.handlers.plotly_buttonclicked({button: {args: [key, value]}});
  }
  await flush();
  assert.equal(a.data[0].type, 'surface');
  assert.deepEqual(a.data[0].z, config.metrics[0].cells[1].linear);
  const single = widget();
  draw(single, {}, {type: 'heatmap', models: 'Model', metrics: [{
    name: 'Estimate', unit: '', years: 2000, ages: 1,
    linear: [-.01, .01], log: [0, 1],
    cells: [{linear: [0], log: [null], text: ['zero']}]
  }]});
  await flush();
  assert.deepEqual(single.data[0].z, [[0]]);
  assert.deepEqual(single.data[0].x, [2000]);
  assert.deepEqual(single.layout.sliders, []);
  single.handlers.plotly_buttonclicked({button: {args: ['scale', 'log']}});
  await flush();
  assert.equal(single.data[0].visible, false);
  console.log('All matrix control combinations passed');
})().catch(error => { console.error(error); process.exit(1); });
