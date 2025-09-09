---
title: Traits
permalink: traits/
layout: base
---

{: .text-center }
### **{{ site.data.stats.n_traits }}** traits &middot; **{{ site.data.stats.n_loci }}** associated loci &middot; **{{ site.data.stats.n_gene_trait_assocs }}**  gene/trait associations

<div id="tag-filter-bar" class="mb-3">
  <div class="small text-muted mb-1">
    <i class="fas fa-tags"></i> Filter by tag
  </div>
  <div id="tag-buttons" class="d-flex flex-wrap gap-2"></div>
  <div id="active-tag-hint" class="small text-muted mt-1" style="display:none;"></div>  
</div>

<div class="table-responsive">
  <table class="table table-hover">
    <thead>
      <tr>
        <th data-bs-toggle="tooltip" data-bs-placement="top" title="Trait">Trait</th>
        <th data-bs-toggle="tooltip" data-bs-placement="top" title="Number of individuals in the GWAS">N</th>
        <th data-bs-toggle="tooltip" data-bs-placement="top" title="Number of associated genomic loci"># loci</th>
        <th data-bs-toggle="tooltip" data-bs-placement="top" title="Number of independent genes selected in the per-locus joint analyses"># indep. genes</th>
        <th data-bs-toggle="tooltip" data-bs-placement="top" title="Total number of genes with significant associations"># total genes</th>
        <th data-bs-toggle="tooltip" data-bs-placement="top" title="Project that produced the GWAS data">Project</th>
        <th data-bs-toggle="tooltip" data-bs-placement="top" title="Link to download the FUSION TWAS association data">Data</th>
      </tr>
    </thead>
    <tbody>
      {% for trait in site.data.traits %}
      <tr data-tags="{{ trait.TAGS }}">
        <td><a href="{{ site.baseurl }}traits/{{ trait.ID }}">{{ trait.NAME }}</a></td>
        <td>{{ trait.N }}</td>
        <td>{{ trait.NUM_LOCI }}</td>
        <td>{{ trait.NUM_JOINT_GENES }}</td>
        <td>{{ trait.NUM_GENES }}</td>
        <td><a href="{{ site.baseurl }}projects/#{{ trait.PROJECT }}">{{ trait.PROJECT }}</a></td>
        <td><a href="{{ site.baseurl }}data/{{ trait.ID }}.tar.bz2"><i class="far fa-file-archive" aria-hidden="true"></i></a></td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
</div>

<script type="text/javascript" class="init">
    $(document).ready(function () {
        const table = $('table').DataTable({
            lengthChange: false,
            paging: false,
            info: false,
            searching: true,
            scrollX: true,
            language: {
                search: '<i class="fa fa-search fa-2x" aria-hidden="true"></i>'
            },
            layout: {
                topStart: function() {
                    let div = document.createElement('div');
                    div.className = 'small';
                    div.innerHTML = '<i class="fas fa-info-circle"></i> Hover over column headers to see descriptions';
                    return div;
                }
            },
            columnDefs: [
                { className: "dt-left dt-head-left", targets: [0, 5, 6] },
                { className: "dt-right dt-head-right", targets: [1, 2, 3, 4] },
            ],
            order: [[2, 'desc']]
        });

        // Build tag list and counts from row data-tags
        const freq = new Map(); // lower-case key -> count
        const displayMap = new Map(); // lower-case key -> display string
        const keys = new Set();
        $('tbody tr').each(function(){
          const t = ($(this).attr('data-tags') || '').split(',').map(s => s.trim()).filter(Boolean);
          if (t.length === 1 && (t[0].toLowerCase() === 'na' || t[0] === '')) return; // skip missing
          t.forEach(tag => {
            const key = tag.toLowerCase();
            keys.add(key);
            if (!displayMap.has(key)) displayMap.set(key, tag);
            freq.set(key, (freq.get(key) || 0) + 1);
          });
        });

        // Sort tags by count desc, then alpha by display label
        const sortedKeys = Array.from(keys).sort((a,b) => {
          const da = -(freq.get(a) || 0), db = -(freq.get(b) || 0);
          if (da !== db) return da - db;
          return (displayMap.get(a) || a).localeCompare(displayMap.get(b) || b);
        });

        // Render buttons
        const $btns = $('#tag-buttons');
        const allBtn = $('<button>', {
          type: 'button',
          class: 'btn btn-sm btn-outline-secondary',
          text: 'All',
          'aria-pressed': 'true'
        });
        $btns.append(allBtn);
        sortedKeys.forEach(key => {
          const tag = displayMap.get(key) || key;
          const count = freq.get(key) || 0;
          const label = `${tag} (${count})`;
          const $b = $('<button>', {
            type: 'button',
            class: 'btn btn-sm btn-outline-primary',
            text: label,
          }).data('tag', tag);
          $btns.append($b);
        });

        let selectedTag = null; // store original-case tag label
        const tableNode = table.table().node();

        // Custom DataTables search for tag filtering (only affects this table)
        $.fn.dataTable.ext.search.push(function(settings, data, dataIndex) {
          if (settings.nTable !== tableNode) return true; // other tables unaffected
          if (!selectedTag) return true;
          const node = table.row(dataIndex).node();
          const rowTags = ($(node).attr('data-tags') || '')
            .split(',')
            .map(s => s.trim().toLowerCase())
            .filter(Boolean);
          return rowTags.includes(selectedTag.toLowerCase());
        });

        const $hint = $('#active-tag-hint');

        function updateHint() {
          if (!selectedTag) {
            $hint.hide().text('');
          } else {
            $hint.show().text(`Showing traits tagged: ${selectedTag}`);
          }
        }

        // Wire up clicks
        $btns.on('click', 'button', function(){
          const tag = $(this).data('tag') || null; // null means All
          if (tag === selectedTag) {
            // toggle off
            selectedTag = null;
          } else {
            selectedTag = tag;
          }
          // Update button active styles
          $btns.find('button').removeClass('active');
          if (selectedTag) {
            $btns.find('button').filter(function(){ return $(this).data('tag') === selectedTag; }).addClass('active');
          } else {
            allBtn.addClass('active');
          }
          updateHint();
          table.draw();
        });

        // Initialize styles
        allBtn.addClass('active');
        updateHint();
    });
</script>
