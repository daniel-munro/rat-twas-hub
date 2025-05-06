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
        <th>Trait</th>
        <th>N</th>
        <th># loci</th>
        <th># indep. genes</th>
        <th># total genes</th>
        <th>Project</th>
        <th>Data</th>
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
        <td><a href="{{ site.baseurl }}projects/">{{ trait.PROJECT }}</a></td>
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
                topStart: 'search',
                topEnd: null,
            },
            columnDefs: [
                { className: "dt-left dt-head-left", targets: [0, 5, 6] },
                { className: "dt-right dt-head-right", targets: [1, 2, 3, 4] },
            ],
            order: [[2, 'desc']]
        });
    });
</script>
