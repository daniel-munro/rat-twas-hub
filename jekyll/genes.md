---
title: Genes
permalink: genes/
layout: default
---

{: .text-center}
### **{{ site.data.genes_stats.n_genes }}** genes &middot; **{{ site.data.genes_stats.n_models }}** models

To search for human genes, try the [cross-species search](#cross-species).

{: #genes}
| Gene | ID | # associated traits | # models |
| --- | --- | --- | --- |
| | | | |

{: #cross-species}
## Cross-species gene search

Each row is an ortholog pair or a gene that has no ortholog in the other species. Rat gene symbols/IDs with links point to the gene entries in this Rat TWAS Hub, while human gene symbols/IDs with links point to the gene entries in the human [TWAS Hub](http://twas-hub.org). Symbols/IDs without links are not included in the TWAS hub of that species.

{: #cross-species .display}
| Human Symbol | Human ID | Rat Symbol | Rat ID | Human Assoc. | Rat Assoc. |
|---|---|---|---|---|---|
| | | | | | |

<script type="text/javascript" class="init">
    $(document).ready(function () {
        $('table#genes').DataTable({
            "ajax": '{{ site.baseurl }}genes.json',
            lengthChange: false,
            "language": {
                search: '<i class="fa fa-search fa-2x" aria-hidden="true"></i>'
            },
            columnDefs: [
                { "targets": [0, 1], "className": "dt-left dt-head-left" },
                { "targets": [2, 3], "className": "dt-right dt-head-right" },
                { "targets": [2, 3], "searchable": false },
            ],
            order: [[2, 'desc']],
        });
        $('table#cross-species').DataTable({
            "ajax": '{{ site.baseurl }}cross-species.json',
            lengthChange: false,
            "language": {
                search: '<i class="fa fa-search fa-2x" aria-hidden="true"></i>'
            },
            columnDefs: [ 
                { "targets": [0, 1, 2, 3], "className": "dt-left dt-head-left" },
                { "targets": [4, 5], "className": "dt-right dt-head-right" },
                { "targets": [4, 5] , "searchable": false },
            ],
            order: [[5, 'desc']],
        });
    });
</script>
