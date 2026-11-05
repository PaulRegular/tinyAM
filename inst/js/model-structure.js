(function() {
  'use strict';
  var root = document.currentScript.closest('.tam-model-structure');
  if (!root) return;
  var scheduled = false;
  var ns = 'http://www.w3.org/2000/svg';
  function coordinate(x) { return Math.round(x * 1000) / 1000; }
  function same(a, b) { return a[0] === b[0] && a[1] === b[1]; }

  function box(element, origin, padding) {
    var r = element.getBoundingClientRect();
    return {left: coordinate(r.left - origin.left - padding),
      right: coordinate(r.right - origin.left + padding),
      top: coordinate(r.top - origin.top - padding),
      bottom: coordinate(r.bottom - origin.top + padding)};
  }
  function port(rect, side, fraction, distance, bounds) {
    var p = (side === 'left' || side === 'right' ?
      [rect[side], rect.top + fraction * (rect.bottom - rect.top)] :
      [rect.left + fraction * (rect.right - rect.left), rect[side]]).map(coordinate);
    var direction = {left: [-1, 0], right: [1, 0], top: [0, -1], bottom: [0, 1]}[side];
    // Keep outer lanes inside the canvas even at large browser font sizes.
    var outside = [Math.max(2, Math.min(bounds.width - 2, p[0] + distance * direction[0])),
      Math.max(2, Math.min(bounds.height - 2, p[1] + distance * direction[1]))].map(coordinate);
    return {tip: p, outside: outside, axis: direction[0] ? 0 : 1};
  }
  function clear(a, b, obstacles) {
    return !obstacles.some(function(r) {
      return a[0] === b[0] ? a[0] > r.left && a[0] < r.right &&
        Math.max(a[1], b[1]) > r.top && Math.min(a[1], b[1]) < r.bottom :
        a[1] > r.top && a[1] < r.bottom && Math.max(a[0], b[0]) > r.left &&
        Math.min(a[0], b[0]) < r.right;
    });
  }
  // Each search state includes its incoming axis. A bend costs 18 CSS pixels,
  // discouraging equal-length staircases while retaining obstacle avoidance.
  function route(from, to, obstacles, width, height) {
    var start = from.outside, end = to.outside;
    if (same(start, end)) return [start];
    function coordinates(values, maximum) {
      return Array.from(new Set(values.map(coordinate))).filter(function(x) {
        return x >= 2 && x <= maximum - 2;
      }).sort(function(a, b) { return a - b; });
    }
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
    var first = id(start) * 2 + from.axis, last = id(end);
    var costs = Array(xs.length * ny * 2).fill(Infinity), turns = costs.slice();
    var previous = Array(costs.length).fill(-1);
    var open = [{id: first, cost: 0, turns: 0, score: distance(start, end)}];
    costs[first] = 0; turns[first] = 0;
    while (open.length) {
      open.sort(function(a, b) {
        return b.score - a.score || b.turns - a.turns || b.cost - a.cost || b.id - a.id;
      });
      var current = open.pop();
      if (Math.abs(current.cost - costs[current.id]) > 1e-7 || current.turns !== turns[current.id]) continue;
      var currentPoint = Math.floor(current.id / 2), axis = current.id % 2;
      if (currentPoint === last) {
        var result = [];
        for (var k = current.id; k !== -1; k = previous[k]) result.push(point(Math.floor(k / 2)));
        return result.reverse();
      }
      var ix = Math.floor(currentPoint / ny), iy = currentPoint % ny;
      [[ix - 1, iy, 0], [ix + 1, iy, 0], [ix, iy - 1, 1], [ix, iy + 1, 1]].forEach(function(next) {
        if (next[0] < 0 || next[0] >= xs.length || next[1] < 0 || next[1] >= ny) return;
        var nextPoint = next[0] * ny + next[1], nextId = nextPoint * 2 + next[2];
        var a = point(currentPoint), b = point(nextPoint);
        if (!clear(a, b, obstacles)) return;
        var bends = Number(axis !== next[2]) + Number(nextPoint === last && next[2] !== to.axis);
        var cost = current.cost + distance(a, b) + 18 * bends;
        var count = current.turns + bends;
        if (cost > costs[nextId] - 1e-7 && !(Math.abs(cost - costs[nextId]) < 1e-7 && count < turns[nextId])) return;
        costs[nextId] = cost; turns[nextId] = count; previous[nextId] = current.id;
        open.push({id: nextId, cost: cost, turns: count, score: cost + distance(b, end)});
      });
    }
    return [];
  }
  function arrow(svg, points, name, dashed) {
    var simple = [];
    points.forEach(function(point) {
      var p = point.map(coordinate);
      if (simple.length && same(simple[simple.length - 1], p)) return;
      while (simple.length > 1) {
        var a = simple[simple.length - 2], b = simple[simple.length - 1];
        if (!(a[0] === b[0] && b[0] === p[0]) && !(a[1] === b[1] && b[1] === p[1])) break;
        simple.pop();
      }
      if (!simple.length || !same(simple[simple.length - 1], p)) simple.push(p);
    });
    if (simple.length < 2) return;
    var group = document.createElementNS(ns, 'g');
    group.setAttribute('data-connection', name);
    var line = document.createElementNS(ns, 'polyline');
    line.setAttribute('points', simple.map(function(p) { return p.join(','); }).join(' '));
    line.setAttribute('fill', 'none'); line.setAttribute('stroke', dashed ? '#af8a27' : '#526a75');
    line.setAttribute('stroke-width', '1.6');
    line.setAttribute('stroke-linejoin', 'round'); line.setAttribute('stroke-linecap', 'round');
    if (dashed) line.setAttribute('stroke-dasharray', '4,3');
    group.appendChild(line);
    var end = simple[simple.length - 1], before = simple[simple.length - 2];
    var length = Math.hypot(end[0] - before[0], end[1] - before[1]);
    var dx = (end[0] - before[0]) / length, dy = (end[1] - before[1]) / length;
    var head = document.createElementNS(ns, 'polygon');
    head.setAttribute('points', [end, [end[0] - 7 * dx - 3 * dy, end[1] - 7 * dy + 3 * dx],
      [end[0] - 7 * dx + 3 * dy, end[1] - 7 * dy - 3 * dx]].map(function(p) { return p.join(','); }).join(' '));
    head.setAttribute('fill', dashed ? '#af8a27' : '#526a75');
    group.appendChild(head); svg.appendChild(group);
  }
  function edgesFor(layout, feedback) {
    var edges;
    if (layout === 'wide') {
      edges = [
        ['initial', 'population', 'top', 'top', .35, .35],
        ['recruitment', 'population', 'right', 'left', .55, .5],
        ['fishing', 'population', 'right', 'bottom', .2, .2],
        ['natural', 'population', 'right', 'bottom', .2, .65],
        ['population', 'survey', 'right', 'left', .4, .5],
        ['population', 'catch', 'right', 'left', .8, .3],
        ['fishing', 'catch', 'bottom', 'bottom', .5, .5],
        ['population', 'ssb', 'bottom', 'top', .5, .5]
      ];
      if (feedback) edges.push(['ssb', 'recruitment', 'left', 'bottom', .5, .8]);
    } else if (layout === 'compact') {
      edges = [
        ['initial', 'population', 'left', 'left', .5, .2, 18],
        ['recruitment', 'population', 'left', 'top', .7, .75],
        ['fishing', 'population', 'left', 'left', .7, .45],
        ['natural', 'population', 'left', 'right', .7, .25],
        ['population', 'survey', 'left', 'left', .7, .35],
        ['population', 'catch', 'right', 'left', .7, .3],
        ['fishing', 'catch', 'left', 'bottom', .85, .75, 26],
        ['population', 'ssb', 'right', 'left', .4, .5]
      ];
      if (feedback) edges.push(['ssb', 'recruitment', 'right', 'right', .4, .8, 26]);
    } else {
      edges = [
        ['initial', 'population', 'left', 'left', .2, .2, 18],
        ['recruitment', 'population', 'right', 'right', .2, .4, 18],
        ['fishing', 'population', 'left', 'left', .65, .6],
        ['natural', 'population', 'right', 'right', .65, .8],
        ['population', 'survey', 'right', 'right', .55, .3],
        ['population', 'catch', 'left', 'left', .75, .35],
        ['fishing', 'catch', 'right', 'right', .85, .75, 26],
        ['population', 'ssb', 'left', 'left', .95, .5, 26]
      ];
      if (feedback) edges.push(['ssb', 'recruitment', 'left', 'left', .8, .8, 26]);
    }
    return edges;
  }
  function draw(graph) {
    var bounds = graph.getBoundingClientRect();
    if (!bounds.width || !bounds.height) return;
    var layout = bounds.width < 560 ? 'stacked' : bounds.width < 1000 ? 'compact' : 'wide';
    graph.setAttribute('data-layout', layout);
    bounds = graph.getBoundingClientRect();
    var extent = {width: coordinate(bounds.width), height: coordinate(bounds.height)};
    var svg = graph.querySelector('.tam-model-connectors'), nodes = {}, obstacles = [];
    graph.querySelectorAll('[data-component]').forEach(function(node) {
      nodes[node.getAttribute('data-component')] = box(node, bounds, 0);
      obstacles.push(box(node, bounds, 4));
    });
    graph.querySelectorAll('[data-input-component]').forEach(function(input) {
      obstacles.push(box(input, bounds, 4));
    });
    svg.replaceChildren();
    svg.setAttribute('viewBox', '0 0 ' + extent.width + ' ' + extent.height);
    edgesFor(layout, graph.getAttribute('data-feedback') === 'true').forEach(function(edge) {
      if (!nodes[edge[0]] || !nodes[edge[1]]) return;
      var lane = edge[6] || 10;
      var from = port(nodes[edge[0]], edge[2], edge[4], lane, extent);
      var to = port(nodes[edge[1]], edge[3], edge[5], lane, extent);
      var middle = route(from, to, obstacles, extent.width, extent.height);
      if (middle.length) arrow(svg, [from.tip].concat(middle, [to.tip]), edge[0] + '-' + edge[1], false);
    });
    graph.querySelectorAll('[data-input-component]').forEach(function(input) {
      var target = nodes[input.getAttribute('data-input-component')];
      if (!target) return;
      var r = box(input, bounds, 0), x = coordinate((target.left + target.right) / 2);
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
