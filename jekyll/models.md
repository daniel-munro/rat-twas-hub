---
title: Models
permalink: models/
---

# Models 

<div class="table-responsive">
  <table class="table table-hover">
    <thead>
      <tr>
        <th>Tissue</th>
        <th>Modality</th>
        <th>Sample size</th>
      </tr>
    </thead>
    <tbody>
      {% for model in site.data.panels %}
      <tr>
        <td>{{ model.TISSUE }}</td>
        <td>{{ model.MODALITY }}</td>
        <td>{{ model.N }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
</div>
