---
title: Projects
permalink: projects/
---

# Projects

Each project provided GWAS data for one or more traits.
Only traits from the original studies with sufficient information and biological interpretability are included in this portal.

<div class="table-responsive">
  <table id="projects" class="table table-hover">
    <thead>
      <tr>
        <th>ID</th>
        <th>Principal investigator</th>
        <th>Title</th>
        <th>Publication</th>
        <th>Animal source</th>
      </tr>
    </thead>
    <tbody>
      {% for project in site.data.projects %}
      <tr id="{{ project.id }}">
        <td>{{ project.id | default: '-' }}</td>
        <td>{{ project.pi | default: '-' }}</td>
        <td>{{ project.title | default: '-' }}</td>
        {% if project.pub_title and project.pub_title != '' %}
        <td><a href="{{ project.pub_url }}">{{ project.pub_title }}</a></td>
        {% else %}
        <td>-</td>
        {% endif %}
        <td>{{ project.animal_source | default: '-' }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
</div>
