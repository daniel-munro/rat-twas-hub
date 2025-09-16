---
title: Models
permalink: models/
---

# Models

A transcriptomic model was trained for each RNA phenotype in each tissue. The table below includes models with cis-heritability P-value &lt; 0.01, which were used for TWAS association testing. Thus, for single phenotype modalities such as expression level, there are 0 or 1 models per gene per tissue, while for multi-phenotype modalities such as isoform ratio, there may be multiple models per gene per tissue.

<div class="table-responsive">
  <table class="table table-hover">
    <thead>
      <tr>
        <th>Tissue</th>
        <th>Modality</th>
        <th>Sample size</th>
        <th>Models</th>
      </tr>
    </thead>
    <tbody>
      {% for model in site.data.panels %}
      <tr>
        <td>{{ model.TISSUE }}</td>
        <td>{{ model.MODALITY }}</td>
        <td>{{ model.N }}</td>
        <td>{{ model.N_MODELS }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
</div>
