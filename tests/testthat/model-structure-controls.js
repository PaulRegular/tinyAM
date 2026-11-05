const fs = require('fs'), vm = require('vm'), assert = require('assert');
const source = fs.readFileSync(process.argv[2], 'utf8');
const order = ['initial', 'recruitment', 'fishing', 'natural', 'population', 'survey', 'catch', 'ssb'];
function element() {
  return {attributes: {}, children: [],
    setAttribute(name, value) { this.attributes[name] = String(value); },
    getAttribute(name) { return this.attributes[name]; },
    appendChild(child) { this.children.push(child); },
    replaceChildren() { this.children = []; }};
}
function widget(width, feedback, options = {}) {
  const queue = [], observations = [], ancestors = [{parentElement: null}];
  const root = {parentElement: ancestors[0], visible: true};
  const pane = {style: {display: ''}, getAttribute: () => '1'};
  const reference = {style: {display: 'none'}, getAttribute: () => '2'};
  const selector = {value: '1', addEventListener: (name, fun) => { selector[name] = fun; }};
  const svg = element(), graph = element(); graph.width = width;
  graph.attributes['data-feedback'] = String(feedback);
  const inputNames = options.inputs || ['survey', 'fishing'];
  const origin = options.fractional ? {left: 73.8125, top: 15.46875} : {left: 0, top: 0};
  const margin = options.fractional ? 28.1372 : 28, gap = options.fractional ? 24.3375 : 24;
  graph.layout = () => {
    const rows = graph.width >= 1000 ?
      [['initial', 'recruitment', 'population', 'survey'], ['fishing', 'natural', 'ssb', 'catch']] :
      graph.width >= 560 ? [['initial', 'recruitment'], ['fishing', 'natural'], ['population', 'ssb'], ['survey', 'catch']] :
      order.map(name => [name]);
    const result = {}, columns = rows[0].length;
    const cellWidth = (graph.width - 2 * margin - (columns - 1) * gap) / columns;
    let top = margin;
    rows.forEach(row => {
      let height = 0;
      row.forEach((name, column) => {
        const nodeHeight = 120 + (options.fractional ? (order.indexOf(name) % 3) * 17.2734 : 0);
        const inputHeight = inputNames.includes(name) ? 42 : 0;
        result[name] = [margin + column * (cellWidth + gap), top + inputHeight, cellWidth, nodeHeight];
        height = Math.max(height, nodeHeight + inputHeight);
      });
      top += height + gap;
    });
    return result;
  };
  function bounds(name, input) {
    const p = graph.layout()[name], top = p[1] - (input ? 42 : 0), height = input ? 26 : p[3];
    return {left: p[0] + origin.left, right: p[0] + p[2] + origin.left,
      top: top + origin.top, bottom: top + height + origin.top, width: p[2], height};
  }
  const nodes = order.map(name => ({getAttribute: () => name, getBoundingClientRect: () => bounds(name, false)}));
  const inputs = inputNames.map(name => ({getAttribute: () => name, getBoundingClientRect: () => bounds(name, true)}));
  graph.getBoundingClientRect = () => ({...origin,
    width: root.visible && pane.style.display !== 'none' ? graph.width : 0,
    height: Math.max(...Object.values(graph.layout()).map(p => p[1] + p[3])) + margin});
  graph.querySelector = () => svg;
  graph.querySelectorAll = query => query === '[data-component]' ? nodes : inputs;
  root.querySelector = () => selector;
  root.querySelectorAll = query => query === '.tam-model-graph' ? [graph] : [pane, reference];
  const observers = {};
  const window = {
    requestAnimationFrame: fun => queue.push(fun), addEventListener: (name, fun) => { window[name] = fun; },
    ResizeObserver: class {constructor(fun) { observers.resize = fun; } observe(target) { observations.push(target); }},
    MutationObserver: class {constructor(fun) { observers.visibility = fun; } observe() {}}
  };
  const document = {currentScript: {closest: () => root}, createElementNS: element};
  vm.runInNewContext(source, {window, document});
  function flush() { while (queue.length) queue.shift()(); }
  function validate() {
    assert.equal(svg.children.length, 8 + inputNames.length + Number(feedback));
    const obstacles = nodes.concat(inputs).map(node => {
      const r = node.getBoundingClientRect();
      return {left: r.left - origin.left, right: r.right - origin.left,
        top: r.top - origin.top, bottom: r.bottom - origin.top};
    });
    svg.children.filter(group => !group.attributes['data-connection'].startsWith('input-')).forEach(group => {
      const points = group.children[0].attributes.points.split(' ').map(p => p.split(',').map(Number));
      assert.equal(group.children[0].attributes['stroke-linejoin'], 'round');
      assert.equal(group.children[0].attributes['stroke-linecap'], 'round');
      assert(points.length <= (graph.width < 560 ? 4 : 8),
        'connector should avoid unnecessary bends: ' + group.attributes['data-connection']);
      points.slice(1).forEach((p, i) => {
        const a = points[i];
        assert(a[0] === p[0] || a[1] === p[1], 'connector must be orthogonal');
        assert(Math.hypot(a[0] - p[0], a[1] - p[1]) > .002, 'connector must not have zero-length segments');
        obstacles.forEach(r => {
          const epsilon = .002; // Subpixel border rounding is not a text crossing.
          const crosses = a[0] === p[0] ? a[0] > r.left + epsilon && a[0] < r.right - epsilon &&
            Math.max(a[1], p[1]) > r.top + epsilon && Math.min(a[1], p[1]) < r.bottom - epsilon :
            a[1] > r.top + epsilon && a[1] < r.bottom - epsilon && Math.max(a[0], p[0]) > r.left + epsilon &&
            Math.min(a[0], p[0]) < r.right - epsilon;
          assert(!crosses, 'connector must avoid text boxes: ' + group.attributes['data-connection']);
        });
      });
    });
    assert(svg.children.some(group => group.children[0].attributes['stroke-dasharray'] === '4,3'));
  }
  flush(); validate();
  assert(observations.includes(root) && observations.includes(graph));
  return {root, graph, svg, selector, pane, reference, observers, flush, validate};
}
const a = widget(1000, true), b = widget(1000, false);
a.graph.width = 640; a.observers.resize(); a.flush(); a.validate();
assert.equal(a.graph.attributes['data-layout'], 'compact');
assert.equal(b.graph.attributes['data-layout'], 'wide');
a.graph.width = 400; a.observers.resize(); a.flush(); a.validate();
assert.equal(a.graph.attributes['data-layout'], 'stacked');
assert(a.svg.attributes.viewBox.startsWith('0 0 400 '));
a.selector.value = '2'; a.selector.change(); a.flush();
assert.equal(a.pane.style.display, 'none'); assert.equal(a.reference.style.display, '');
assert.equal(b.pane.style.display, '');
a.graph.width = 1000; a.selector.value = '1'; a.selector.change(); a.flush(); a.validate();
assert.equal(a.graph.attributes['data-layout'], 'wide');
a.root.visible = false; a.graph.width = 640; a.observers.visibility(); a.flush();
a.root.visible = true; a.observers.visibility(); a.flush(); a.validate();
assert.equal(a.graph.attributes['data-layout'], 'compact');
assert(a.graph.getBoundingClientRect().height < 900, 'compact layout should use four filled rows');

const wide = widget(1200, true);
assert(wide.graph.getBoundingClientRect().height < 500, 'desktop layout should use two filled rows');
const ssb = wide.svg.children.find(group => group.attributes['data-connection'] === 'population-ssb');
assert.equal(ssb.children[0].attributes.points.split(' ').length, 2, 'aligned population/SSB arrow should be straight');

// Fractional browser pixels and uneven cards exposed diagonal subpixel stubs.
const fractional = widget(1100.675, true, {fractional: true,
  inputs: ['recruitment', 'population', 'fishing', 'natural', 'catch', 'survey']});
fractional.graph.width = 760.287; fractional.observers.resize(); fractional.flush(); fractional.validate();
fractional.graph.width = 430.391; fractional.observers.resize(); fractional.flush(); fractional.validate();
console.log('Responsive model structure controls passed');
