---
title: Cross-species gene search
permalink: cross-species/
layout: default
---

# Cross-species gene search

Search and compare gene associations across species:

| Human Symbol | Human ID | Rat Symbol | Rat ID | Human Assoc. | Rat Assoc. |
|---|---|---|---|---|---|
| | | | | | |
{: #genes .display}

<script type="text/javascript" class="init">
    $(document).ready(function() {
        $('table#genes').DataTable( {
            "ajax": '{{ site.baseurl }}cross-species.json',
            lengthChange: false,
            "language": {
                search: '<i class="fa fa-search fa-2x" aria-hidden="true"></i>'
            },
            columnDefs: [ { "targets": [4,5] , "searchable": false }, { "className": "dt-left dt-head-left", "targets": [0,1,2,3] },{ "className": "dt-right", "targets": '_all' } ],
            order: [[5, 'desc']],
        } );	  						
    } );
</script>
