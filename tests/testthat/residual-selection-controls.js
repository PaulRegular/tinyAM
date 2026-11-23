const fs = require('fs'), assert = require('assert');
const initialize = eval('(' + fs.readFileSync(process.argv[2], 'utf8') + ')');
const configs = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
const config = configs.osa;

class Element {
  constructor(tag) {
    this.tagName = tag;
    this.children = [];
    this.handlers = {};
    this.style = {};
  }
  appendChild(child) { this.children.push(child); }
  replaceChildren() { this.children = []; }
  get childElementCount() { return this.children.length; }
  addEventListener(name, listener) { this.handlers[name] = listener; }
  choose(value) {
    assert(this.children.some(option => option.value === String(value)));
    this.value = String(value);
    this.handlers.change();
  }
  find(className) {
    if (this.className === className) return this;
    return this.children.map(child => child.find(className)).find(Boolean);
  }
}
global.document = {createElement: name => new Element(name)};
global.window = {getComputedStyle: () => ({position: 'static'})};
global.Plotly = {
  react: (el, data, layout) => new Promise(resolve => {
    assert.equal(el.draw, undefined, 'Widget runs overlapping Plotly updates');
    el.draw = () => {
      el.data = data;
      el.layout = layout;
      el.draw = undefined;
      resolve();
    };
    el.nDraws = (el.nDraws || 0) + 1;
  })
};
const flush = () => new Promise(resolve => setImmediate(resolve));
const finish = async el => { await flush(); el.draw(); await flush(); };
const points = el => el.data.reduce((n, trace) => n + trace.x.length, 0);
const years = el => el.data.flatMap(trace => trace.x).sort();
function widget(options = config) {
  const el = new Element('div');
  initialize(el, {layout: {}}, JSON.parse(JSON.stringify(options)));
  return el;
}

async function verifyObservations() {
  const a = widget(), b = widget();
  await finish(a);
  await finish(b);
  assert.equal(points(a), 2);
  assert.equal(points(b), 2);
  const residual = a.find('tam-residual-type'), model = a.find('tam-residual-frame');
  assert.deepEqual(residual.children.map(option => option.textContent), ['Standardized / Quantile', 'OSA']);
  assert.deepEqual(model.children.map(option => option.textContent), ['A', 'B']);
  assert.equal(residual.value, '1');
  assert.equal(model.value, '0');
  assert.equal(a.layout.margin.t, 60);
  assert.equal(a.layout.title, '');

  model.choose(1);
  assert.equal(residual.value, '0'); // Auto uses the next model's default.
  await finish(a);
  assert.equal(points(a), 6);
  residual.choose(1);
  await finish(a);
  assert.equal(points(a), 0);
  assert(a.layout.annotations[0].text.includes('OSA residuals unavailable'));
  assert.equal(points(b), 2);
  assert.equal(b.find('tam-residual-frame').value, '0');

  model.choose(0);
  assert.equal(residual.value, '1'); // Explicit selection persists.
  await finish(a);
  assert.equal(points(a), 2);
  assert.deepEqual(a.layout.annotations, []);
  assert(a.data.some(trace => trace.text.some(label => label.includes('oneStepGeneric'))));
  assert.deepEqual(years(a), [2000, 2001]);

  // Keep one render artificially pending while several selections arrive.
  residual.choose(0);
  await flush();
  const before = a.nDraws;
  model.choose(1);
  residual.choose(1);
  model.choose(0);
  model.choose(1);
  await flush();
  assert.equal(a.nDraws, before);
  await finish(a);
  assert.equal(a.nDraws, before + 1);
  await finish(a);
  assert.equal(points(a), 0);
  residual.choose(0);
  await finish(a);
  assert.equal(points(a), 6);
  assert.equal(a.find('tam-residual-type'), residual);
  assert.equal(a.find('tam-residual-frame'), model);
  assert.deepEqual(years(a), [2000, 2000, 2001, 2001, 2002, 2002]);
  assert.equal(points(b), 2);

  const one = JSON.parse(JSON.stringify(config));
  one.frames = 'Only';
  one.showSelector = false;
  one.views.forEach(view => { view.cells = [view.cells[0]]; });
  one.views[one.view].cells[0].data = [{x: 2000, y: 1, customdata: 'one', marker: {size: 1}}];
  const single = widget(one);
  await finish(single);
  assert.deepEqual(single.data[0].x, [2000]);
  assert.deepEqual(single.data[0].customdata, ['one']);
  assert.equal(single.childElementCount, 0);
  assert.equal(single.layout.margin.t, 35);

  const delta = widget(configs.delta);
  await finish(delta);
  const deltaType = delta.find('tam-residual-type');
  const deltaModel = delta.find('tam-residual-frame');
  assert.deepEqual(deltaType.children.map(option => option.textContent), ['Standardized / Quantile']);
  deltaModel.choose(1);
  assert.deepEqual(deltaType.children.map(option => option.textContent), ['Standardized', 'Quantile']);
  assert.deepEqual(deltaType.children.map(option => option.value), ['0', '2']);
  assert.equal(deltaType.value, '2');
  await finish(delta);
  assert.equal(points(delta), 6); // Includes the sampled zero.
  deltaType.choose(2);
  await finish(delta);
  deltaModel.choose(0);
  assert.deepEqual(deltaType.children.map(option => option.textContent), ['Standardized / Quantile']);
  assert.equal(deltaType.value, '0');
  await finish(delta);
  deltaModel.choose(1);
  await finish(delta);
  assert.equal(deltaType.value, '2'); // The combined label preserves explicit Quantile.
  assert.equal(points(delta), 6);
  deltaType.choose(0);
  await finish(delta);
  assert.equal(points(delta), 5);
  deltaModel.choose(0);
  await finish(delta);
  deltaModel.choose(1);
  await finish(delta);
  assert.equal(deltaType.value, '0'); // Standardized remains explicit.
  assert.equal(points(delta), 5);

  const sampled = widget(configs.sampled);
  await finish(sampled);
  const sampledType = sampled.find('tam-residual-type');
  const sampledModel = sampled.find('tam-residual-frame');
  assert.deepEqual(sampledType.children.map(option => option.textContent), ['Standardized / Quantile']);
  sampledModel.choose(1);
  await finish(sampled);
  assert.deepEqual(sampledType.children.map(option => option.textContent), ['Standardized', 'Quantile']);
  assert.equal(sampledType.value, '2');
  assert.equal(points(sampled), 5); // The gap is not filled with a fitted-state residual.
  assert(sampled.data.every(trace => trace.text.every(label => label.includes('mle-mvn'))));
  sampledType.choose(2);
  await finish(sampled);
  sampledModel.choose(0);
  await finish(sampled);
  assert.equal(sampledType.value, '0');
  assert.equal(points(sampled), 6);
  assert(sampled.data.every(trace => trace.text.every(label => label.includes('mle-eb) [fallback]'))));
  sampledModel.choose(1);
  await finish(sampled);
  assert.equal(sampledType.value, '2');
  assert.equal(points(sampled), 5);
  console.log('Residual controls passed');
}

async function verifyProcess(config) {
  const a = widget(config), b = widget(config);
  await finish(a);
  await finish(b);
  const process = a.find('tam-residual-type'), model = a.find('tam-residual-frame');
  assert.deepEqual(process.children.map(option => option.textContent), ['N', 'F', 'M']);
  assert.deepEqual(model.children.map(option => option.textContent), ['scaa', 'sam', 'ncam']);
  assert.equal(process.value, '0');
  assert.equal(model.value, '1');
  assert.equal(points(a), 6);
  process.choose(2);
  await finish(a);
  assert.deepEqual(a.layout.yaxis.ticktext, ['2-4', '5-8']);
  model.choose(2);
  await finish(a);
  assert.equal(process.value, '2');
  assert.deepEqual(a.layout.yaxis.ticktext, ['2']);
  assert.deepEqual(a.layout.yaxis.tickvals, [1]);
  assert.equal(points(a), 3);
  assert(a.data.every(trace => trace.text.every(label => label.includes('Model: ncam'))));
  model.choose(0);
  await finish(a);
  assert.equal(points(a), 0);
  assert.equal(a.layout.annotations[0].text, 'Unavailable for this model component');

  process.choose(1);
  await flush();
  model.choose(1);
  process.choose(0);
  model.choose(2);
  await finish(a);
  await finish(a);
  assert.equal(points(a), 6);
  assert(a.data.every(trace => trace.text.every(label => label.includes('Model: ncam<br>N residual'))));
  assert.deepEqual(a.layout.annotations, []);
  assert.equal(b.find('tam-residual-type').value, '0');
  assert.equal(b.find('tam-residual-frame').value, '1');
  console.log('Process bubble controls passed');
}

(configs.process ? verifyProcess(configs.process) : verifyObservations())
  .catch(error => { console.error(error); process.exit(1); });
