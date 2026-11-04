(function() {
  'use strict';
  var root = document.currentScript.closest('.tam-model-structure');
  if (!root) return;
  var scheduled = false;
  var ns = 'http://www.w3.org/2000/svg';

  function box(element, origin, padding) {
    var r = element.getBoundingClientRect();
    return {left: r.left - origin.left - padding, right: r.right - origin.left + padding,
      top: r.top - origin.top - padding, bottom: r.bottom - origin.top + padding};
  }
  function port(rect, side, fraction) {
    var p = side === 'left' || side === 'right' ?
      [rect[side], rect.top + fraction * (rect.bottom - rect.top)] :
      [rect.left + fraction * (rect.right - rect.left), rect[side]];
    var direction = {left: [-1, 0], right: [1, 0], top: [0, -1], bottom: [0, 1]}[side];
    return {tip: p, outside: [p[0] + 10 * direction[0], p[1] + 10 * direction[1]]};
  }
  function clear(a, b, obstacles) {
    return !obstacles.some(function(r) {
      return a[0] === b[0] ? a[0] > r.left && a[0] < r.right &&
        Math.max(a[1], b[1]) > r.top && Math.min(a[1], b[1]) < r.bottom :
        a[1] > r.top && a[1] < r.bottom && Math.max(a[0], b[0]) > r.left &&
        Math.min(a[0], b[0]) < r.right;
    });
  }
  // Orthogonal paths follow measured gaps, keeping arrows outside node text at
  // every layout width. The small visibility grid needs no layout dependency.
  function route(start, end, obstacles, width, height) {
    function round(x) { return Math.round(x * 1000) / 1000; }
    function coordinates(values, maximum) {
      return Array.from(new Set(values.map(round))).filter(function(x) {
        return x >= 2 && x <= maximum - 2;
      }).sort(function(a, b) { return a - b; });
    }
    start = start.map(round); end = end.map(round);
    var xs = coordinates([2, width - 2, start[0], end[0]].concat(obstacles.flatMap(function(r) {
      return [r.left - 2, r.right + 2];
    })), width);
    var ys = coordinates([2, height - 2, start[1], end[1]].concat(obstacles.flatMap(function(r) {
      return [r.top - 2, r.bottom + 2];
    })), height);
    var ny = ys.length;
    function id(p) { return xs.indexOf(p[0]) * ny + ys.indexOf(p[1]); }
    function point(i) { return [xs[Math.floor(i / ny)], ys[i % ny]]; }
    function distance(a, b) { return Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]); }
    if (!xs.includes(start[0]) || !xs.includes(end[0]) || !ys.includes(start[1]) || !ys.includes(end[1])) return [];
    var first = id(start), last = id(end), costs = Array(xs.length * ny).fill(Infinity);
    var previous = Array(costs.length).fill(-1), open = [{id: first, cost: 0, score: distance(start, end)}];
    costs[first] = 0;
    while (open.length) {
      open.sort(function(a, b) { return b.score - a.score || a.cost - b.cost; });
      var current = open.pop();
      if (current.cost !== costs[current.id]) continue;
      if (current.id === last) {
        var result = [];
        for (var k = last; k !== -1; k = previous[k]) result.push(point(k));
        return result.reverse();
      }
      var ix = Math.floor(current.id / ny), iy = current.id % ny;
      [[ix - 1, iy], [ix + 1, iy], [ix, iy - 1], [ix, iy + 1]].forEach(function(next) {
        if (next[0] < 0 || next[0] >= xs.length || next[1] < 0 || next[1] >= ny) return;
        var nextId = next[0] * ny + next[1], a = point(current.id), b = point(nextId);
        if (!clear(a, b, obstacles)) return;
        var cost = current.cost + distance(a, b);
        if (cost >= costs[nextId]) return;
        costs[nextId] = cost; previous[nextId] = current.id;
        open.push({id: nextId, cost: cost, score: cost + distance(b, end)});
      });
    }
    return [];
  }
  function arrow(svg, points, name, dashed) {
    if (points.length < 2) return;
    points = points.filter(function(p, i) {
      if (i === 0 || i === points.length - 1) return true;
      var a = points[i - 1], b = points[i + 1];
      return !(a[0] === p[0] && p[0] === b[0]) && !(a[1] === p[1] && p[1] === b[1]);
    });
    var group = document.createElementNS(ns, 'g');
    group.setAttribute('data-connection', name);
    var line = document.createElementNS(ns, 'polyline');
    line.setAttribute('points', points.map(function(p) { return p.join(','); }).join(' '));
    line.setAttribute('fill', 'none'); line.setAttribute('stroke', dashed ? '#af8a27' : '#526a75');
    line.setAttribute('stroke-width', '1.6');
    if (dashed) line.setAttribute('stroke-dasharray', '4,3');
    group.appendChild(line);
    var end = points[points.length - 1], before = points[points.length - 2];
    var length = Math.hypot(end[0] - before[0], end[1] - before[1]);
    if (!length) return;
    var dx = (end[0] - before[0]) / length, dy = (end[1] - before[1]) / length;
    var head = document.createElementNS(ns, 'polygon');
    head.setAttribute('points', [end, [end[0] - 7 * dx - 3 * dy, end[1] - 7 * dy + 3 * dx],
      [end[0] - 7 * dx + 3 * dy, end[1] - 7 * dy - 3 * dx]].map(function(p) { return p.join(','); }).join(' '));
    head.setAttribute('fill', dashed ? '#af8a27' : '#526a75');
    group.appendChild(head); svg.appendChild(group);
  }
  function draw(graph) {
    var bounds = graph.getBoundingClientRect();
    if (!bounds.width || !bounds.height) return;
    graph.setAttribute('data-layout', bounds.width < 460 ? 'stacked' : bounds.width < 760 ? 'compact' : 'wide');
    bounds = graph.getBoundingClientRect();
    var svg = graph.querySelector('.tam-model-connectors'), nodes = {}, obstacles = [];
    graph.querySelectorAll('[data-component]').forEach(function(node) {
      nodes[node.getAttribute('data-component')] = box(node, bounds, 0);
      obstacles.push(box(node, bounds, 4));
    });
    graph.querySelectorAll('[data-input-component]').forEach(function(input) {
      obstacles.push(box(input, bounds, 4));
    });
    svg.replaceChildren();
    svg.setAttribute('viewBox', '0 0 ' + bounds.width + ' ' + bounds.height);
    var edges = [
      ['initial', 'population', 'right', 'left', .5, .2],
      ['recruitment', 'population', 'right', 'left', .5, .5],
      ['fishing', 'population', 'top', 'bottom', .5, .25],
      ['natural', 'population', 'top', 'left', .5, .8],
      ['population', 'survey', 'right', 'left', .3, .5],
      ['population', 'catch', 'right', 'left', .7, .5],
      ['fishing', 'catch', 'bottom', 'bottom', .5, .5],
      ['population', 'ssb', 'bottom', 'top', .7, .5]
    ];
    if (graph.getAttribute('data-feedback') === 'true') edges.push(['ssb', 'recruitment', 'left', 'right', .5, .8]);
    edges.forEach(function(edge) {
      if (!nodes[edge[0]] || !nodes[edge[1]]) return;
      var from = port(nodes[edge[0]], edge[2], edge[4]), to = port(nodes[edge[1]], edge[3], edge[5]);
      var middle = route(from.outside, to.outside, obstacles, bounds.width, bounds.height);
      if (middle.length) arrow(svg, [from.tip].concat(middle, [to.tip]), edge[0] + '-' + edge[1], false);
    });
    graph.querySelectorAll('[data-input-component]').forEach(function(input) {
      var target = nodes[input.getAttribute('data-input-component')];
      if (!target) return;
      var r = box(input, bounds, 0), x = (target.left + target.right) / 2;
      arrow(svg, [[x, r.bottom], [x, target.top]], 'input-' + input.getAttribute('data-input-component'), true);
    });
  }
  function schedule() {
    if (scheduled) return;
    scheduled = true;
    window.requestAnimationFrame(function() {
      scheduled = false;
      root.querySelectorAll('.tam-model-graph').forEach(draw);
    });
  }
  var selector = root.querySelector('[data-model-selector]');
  if (selector) selector.addEventListener('change', function() {
    var value = selector.value;
    root.querySelectorAll('[data-model-pane]').forEach(function(pane) {
      pane.style.display = pane.getAttribute('data-model-pane') === value ? '' : 'none';
    });
    schedule();
  });
  if (window.ResizeObserver) {
    var resize = new window.ResizeObserver(schedule);
    resize.observe(root);
    root.querySelectorAll('.tam-model-graph').forEach(function(graph) { resize.observe(graph); });
  }
  // Dashboard tabs change ancestor classes/styles without a window resize.
  if (window.MutationObserver) {
    var visibility = new window.MutationObserver(schedule);
    for (var parent = root.parentElement; parent; parent = parent.parentElement) {
      visibility.observe(parent, {attributes: true, attributeFilter: ['class', 'style', 'hidden']});
    }
  }
  window.addEventListener('resize', schedule);
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(schedule);
  schedule();
})();
