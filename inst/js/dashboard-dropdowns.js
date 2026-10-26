(function() {
  function install() {
    document.querySelectorAll('[data-dropdown-tabs]').forEach(function(row) {
      var nav = row.querySelector(':scope > .nav-tabs');
      if (!nav || nav.dataset.dropdownReady) return;
      nav.dataset.dropdownReady = 'true';
      var groups = {};
      row.querySelectorAll(':scope > .tab-content > [data-tab-group]').forEach(function(pane) {
        var name = pane.dataset.tabGroup;
        var anchor = Array.from(nav.querySelectorAll('a[role="tab"]')).find(function(a) {
          return a.getAttribute('href') === '#' + pane.id;
        });
        if (!anchor) return;
        var group = groups[name];
        if (!group) {
          var li = document.createElement('li');
          li.className = 'nav-item dropdown';
          // Match flexdashboard's tab headers so its Bootstrap rules provide
          // the same font, padding and border.
          var toggle = document.createElement('a');
          toggle.className = 'nav-link dropdown-toggle';
          toggle.setAttribute('href', '#');
          toggle.setAttribute('role', 'button');
          toggle.textContent = name;
          toggle.setAttribute('data-toggle', 'dropdown');
          toggle.setAttribute('aria-expanded', 'false');
          toggle.setAttribute('aria-haspopup', 'true');
          toggle.addEventListener('click', function(event) { event.preventDefault(); });
          var menu = document.createElement('div');
          menu.className = 'dropdown-menu';
          li.appendChild(toggle); li.appendChild(menu);
          nav.insertBefore(li, anchor.parentElement);
          group = groups[name] = {toggle: toggle, menu: menu, items: []};
          menu.addEventListener('keydown', function(event) {
            var items = Array.from(menu.querySelectorAll('button'));
            var index = items.indexOf(document.activeElement);
            if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
              event.preventDefault(); event.stopPropagation();
              items[(index + (event.key === 'ArrowDown' ? 1 : -1) + items.length) % items.length].focus();
            } else if (event.key === 'Escape') {
              event.preventDefault(); event.stopPropagation();
              jQuery(toggle).dropdown('hide'); toggle.focus();
            }
          });
        }
        var item = document.createElement('button');
        item.type = 'button'; item.className = 'dropdown-item';
        item.textContent = pane.dataset.tabLabel || anchor.textContent;
        item.setAttribute('aria-controls', pane.id);
        item.addEventListener('click', function() {
          // Activate the existing flexdashboard tab, retaining its chart container.
          jQuery(anchor).tab('show');
          jQuery(group.toggle).dropdown('hide');
          group.toggle.focus();
        });
        group.menu.appendChild(item);
        group.items.push({item: item, anchor: anchor});
        if (!group.selected) group.selected = item;
        anchor.parentElement.hidden = true;
      });
      function sync() {
        Object.keys(groups).forEach(function(name) {
          var group = groups[name], active = false;
          group.items.forEach(function(entry) {
            if (entry.anchor.classList.contains('active')) group.selected = entry.item;
          });
          group.items.forEach(function(entry) {
            var selected = entry.anchor.classList.contains('active');
            entry.item.classList.toggle('dropdown-selected', entry.item === group.selected);
            entry.item.setAttribute('aria-current', selected ? 'true' : 'false');
            active = active || selected;
          });
          // Reserve Bootstrap's active class for its original tab anchors.
          group.toggle.classList.toggle('dropdown-selected', active);
        });
      }
      jQuery(nav).on('shown.bs.tab', function(event) {
        sync();
        var pane = document.getElementById(event.target.getAttribute('aria-controls'));
        requestAnimationFrame(function() {
          if (pane && window.Plotly) pane.querySelectorAll('.js-plotly-plot').forEach(function(plot) {
            if (plot.offsetWidth && plot.offsetHeight) Plotly.Plots.resize(plot);
          });
        });
      });
      sync();
    });
  }
  jQuery(document).on('flexdashboard:layoutcomplete', install);
  jQuery(install);
})();
