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
function widget(width, feedback) {
  const queue = [], observations = [], ancestors = [{parentElement: null}];
  const root = {parentElement: ancestors[0], visible: true};
  const pane = {style: {display: ''}, getAttribute: () => '1'};
  const reference = {style: {display: 'none'}, getAttribute: () => '2'};
  const selector = {value: '1', addEventListener: (name, fun) => { selector[name] = fun; }};
  const svg = element(), graph = element(); graph.width = width;
  graph.attributes['data-feedback'] = String(feedback);
  graph.layout = () => {
    if (graph.width >= 760) return {
      initial: [24, 60, 200, 120], recruitment: [24, 280, 200, 120],
      fishing: [24, 500, 200, 120], natural: [274, 500, 200, 120],
      population: [524, 280, 200, 120], survey: [774, 60, 200, 120],
      catch: [774, 500, 200, 120], ssb: [524, 500, 200, 120]};
    if (graph.width >= 460) return {
      initial: [24, 60, 270, 120], recruitment: [346, 60, 270, 120],
      fishing: [24, 280, 270, 120], natural: [346, 280, 270, 120],
      population: [24, 500, 592, 120], survey: [24, 720, 270, 120],
      catch: [346, 720, 270, 120], ssb: [24, 940, 592, 120]};
    return Object.fromEntries(order.map((name, i) => [name, [24, 60 + i * 200, graph.width - 48, 120]]));
  };
  function bounds(name, input) {
    const p = graph.layout()[name], top = p[1] - (input ? 42 : 0), height = input ? 26 : p[3];
    return {left: p[0], right: p[0] + p[2], top, bottom: top + height, width: p[2], height};
  }
  const nodes = order.map(name => ({getAttribute: () => name, getBoundingClientRect: () => bounds(name, false)}));
  const inputs = ['survey', 'fishing'].map(name => ({getAttribute: () => name, getBoundingClientRect: () => bounds(name, true)}));
  graph.getBoundingClientRect = () => ({left: 0, top: 0,
    width: root.visible && pane.style.display !== 'none' ? graph.width : 0,
    height: graph.width >= 760 ? 660 : graph.width >= 460 ? 1100 : 1640});
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
    assert.equal(svg.children.length, 10 + Number(feedback));
    const obstacles = nodes.concat(inputs).map(node => node.getBoundingClientRect());
    svg.children.filter(group => !group.attributes['data-connection'].startsWith('input-')).forEach(group => {
      const points = group.children[0].attributes.points.split(' ').map(p => p.split(',').map(Number));
      points.slice(1).forEach((p, i) => {
        const a = points[i];
        assert(a[0] === p[0] || a[1] === p[1], 'connector must be orthogonal');
        obstacles.forEach(r => {
          const crosses = a[0] === p[0] ? a[0] > r.left && a[0] < r.right &&
            Math.max(a[1], p[1]) > r.top && Math.min(a[1], p[1]) < r.bottom :
            a[1] > r.top && a[1] < r.bottom && Math.max(a[0], p[0]) > r.left &&
            Math.min(a[0], p[0]) < r.right;
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
console.log('Responsive model structure controls passed');
