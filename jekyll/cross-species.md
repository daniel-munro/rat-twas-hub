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


<!-- <script type="text/javascript">
$(document).ready(function() {
    $('#genes').DataTable({
        "ajax": '{{ site.baseurl }}cross-species.json',
        "columns": [
            { 
                "data": "human_symbol",
                "render": function(data, type, row) {
                    return data || 'N/A';
                }
            },
            { 
                "data": "human_id",
                "render": function(data, type, row) {
                    if (row.human_has_model) {
                        return '<a href="https://portal.brain-map.org/genes/' + data + '">' + data + '</a>';
                    }
                    return data || 'N/A';
                }
            },
            { 
                "data": "rat_symbol",
                "render": function(data, type, row) {
                    return data || 'N/A';
                }
            },
            { 
                "data": "rat_id",
                "render": function(data, type, row) {
                    if (row.rat_has_model) {
                        return '<a href="{{ site.baseurl }}genes/' + data + '">' + data + '</a>';
                    }
                    return data || 'N/A';
                }
            },
            { "data": "human_assoc" },
            { "data": "rat_assoc" }
        ],
        "pageLength": 25,
        "order": [[0, 'asc']],
        "search": {
            "regex": true,
            "smart": true
        }
    });
});
</script>
 -->
