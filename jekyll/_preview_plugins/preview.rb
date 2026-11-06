# frozen_string_literal: true

require "json"

module PortalPreview
  class Reader < Jekyll::Reader
    def initialize(site)
      super
      @traits = read_ids("traits")
      @genes = read_ids("genes")
      @original_include = site.include
      site.exclude += %w[traits/ genes/ assets/images/porcupine/ genes.json cross-species.json]
    end

    def read
      @data_reader = Jekyll::DataReader.new(site)
      # Jekyll reads explicitly included files even inside excluded directories.
      site.include = @original_include + selected_files
      super
      add_gene_index("genes.json", 0, "./")
      add_gene_index("cross-species.json", 2, "../genes/")
    end

    def read_data
      site.data = @preview_data
    end

    private

    def read_ids(key)
      ids = site.config.fetch("preview")[key]
      unless ids.is_a?(Array) && !ids.empty? && ids.all? { |id| id.is_a?(String) && id.match?(/\A[\w-][\w.-]*\z/) }
        raise "Preview #{key} must be a nonempty list of IDs."
      end

      ids.uniq
    end

    def input(relative_path)
      path = site.in_source_dir(relative_path)
      return path if File.file?(path)

      raise "Preview input is missing: #{relative_path}. Check the configured IDs and generate the site inputs first."
    end

    def data(relative_path)
      @data_reader.read_data_file(input("_data/#{relative_path}"))
    end

    def optional_data(relative_path)
      path = site.in_source_dir("_data/#{relative_path}")
      @data_reader.read_data_file(path) if File.file?(path)
    end

    def selected_files
      traits = data("traits.tsv").select { |row| @traits.include?(row["ID"]) }
      missing = @traits - traits.map { |row| row.fetch("ID") }
      raise "Unknown preview trait IDs in _data/traits.tsv: #{missing.join(', ')}" unless missing.empty?

      @preview_data = {
        "traits" => traits,
        "panels" => data("panels.tsv"),
        "projects" => data("projects.tsv"),
        "stats" => data("stats.yml"),
        "gene_models" => {},
        "gene_summaries" => {},
        "trait_hits" => {},
        "trait_loci" => {},
        "trait_panels" => {},
        "trait_pleio" => {},
        "locus_models" => {},
      }

      files = []
      @traits.each do |id|
        files.concat(["traits/#{id}.md", "assets/images/porcupine/#{id}.png"])
        @preview_data["trait_hits"][id] = optional_data("trait_hits/#{id}.tsv")
        @preview_data["trait_panels"][id] = data("trait_panels/#{id}.tsv")
        @preview_data["trait_pleio"][id] = data("trait_pleio/#{id}.tsv")
          .select { |row| @traits.include?(row["trait_id"]) }
        loci = data("trait_loci/#{id}.tsv")
        @preview_data["trait_loci"][id] = loci
        @preview_data["locus_models"][id] = loci.to_h do |locus|
          number = locus.fetch("locus_num")
          files.concat(["traits/#{id}/#{number}.md", "traits/#{id}/#{number}.cond.csv"])
          export = "traits/#{id}/#{number}.models.tsv"
          files << export if File.file?(site.in_source_dir(export))
          [number, data("locus_models/#{id}/#{number}.tsv")]
        end
      end
      @genes.each do |id|
        files << "genes/#{id}.md"
        key = @data_reader.sanitize_filename(id)
        summary = optional_data("gene_summaries/#{id}.json")
        if summary
          summary["hits"].select! { |hit| @traits.include?(hit["trait_id"]) }
          @preview_data["gene_summaries"][key] = summary
        end
        @preview_data["gene_models"][key] = data("gene_models/#{id}.tsv")
      end
      files.each { |file| input(file) }
    end

    def add_gene_index(filename, column, prefix)
      index = JSON.parse(File.read(input(filename)))
      index.fetch("data").select! do |row|
        href = row.fetch(column)[/href="#{Regexp.escape(prefix)}([^"]+)"/, 1]
        @genes.include?(href)
      end
      page = Jekyll::PageWithoutAFile.new(site, site.source, "", filename)
      page.data.merge!("layout" => "none", "render_with_liquid" => false)
      page.content = JSON.generate(index)
      site.pages << page
    end
  end

  Jekyll::Hooks.register :site, :after_init do |site|
    next unless site.config["preview"]

    unless site.dest == site.in_source_dir("_site-preview")
      raise "The portal preview must use _site-preview as its destination."
    end

    site.reader = Reader.new(site)
  end

  Jekyll::Hooks.register :pages, :post_render do |page|
    preview = page.site.config["preview"]
    next unless preview && page.output_ext == ".html"

    # Keep the scientific results visible, but unlink details outside the sample.
    baseurl = Regexp.escape(page.site.baseurl)
    detail_link = %r{<a\b[^>]*\shref\s*=\s*(["'])#{baseurl}(traits|genes)/([^/"'#?]+)(?:/[^"'#?]*)?\1[^>]*>(.*?)</a>}m
    page.output = page.output.gsub(detail_link) do |link|
      preview.fetch(Regexp.last_match(2)).include?(Regexp.last_match(3)) ? link : Regexp.last_match(4)
    end
  end
end
