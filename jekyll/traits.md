---
title: Traits
permalink: traits/
layout: base
---

{: .text-center }
### **{{ site.data.stats.n_traits }}** traits &middot; **{{ site.data.stats.n_loci }}** associated loci &middot; **{{ site.data.stats.n_gene_trait_assocs }}**  gene/trait associations

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
      <tr>
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
        $('table').DataTable({
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
    });
</script>
