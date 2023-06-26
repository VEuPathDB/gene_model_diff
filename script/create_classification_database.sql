drop table if exists exon;
create table if not exists exon(
			exon_id  varchar(500),
			source   varchar(500),
			scaffold varchar(500),
			strand   varchar(500),
			start    int,
			end      int			
);
create index idx_exon1 ON exon(source(1), scaffold(16), strand, start, end);

drop table if exists cds;
create table if not exists cds(
			cds_parent_id varchar(500),
			start         	int,
			end           	int,
			md5_checksum  	varchar(50),
			cds_error_code	int
);
create index idx_cds_parentid on cds(cds_parent_id);

drop table if exists gene_model;
create table if not exists gene_model(
			id 			  int not null AUTO_INCREMENT PRIMARY KEY,
			exon_id       varchar(500),
			transcript_id varchar(500),
			gene_id       varchar(500),
			source        varchar(500),
			error_code	  int,
			biotype       varchar(500)
);
create index idx_gene_model_trid on gene_model(transcript_id, source(1));
create index idx_genemodel_geneid on gene_model(gene_id, source(1));
create index idx_genemodel_exonid on gene_model(exon_id, source(1));

drop table if exists exon_mappings;
create table if not exists exon_mappings(
			id       int not null AUTO_INCREMENT PRIMARY KEY,
			cap_exon_id   varchar(500),
			vb_exon_id    varchar(500),
			map_type      varchar(500)
);
create index idx_exonmap_capid on exon_mappings(cap_exon_id);
create index idx_exonmap_vbid on exon_mappings(vb_exon_id);

drop table if exists gene_clusters;
create table if not exists gene_clusters(
			gene_cluster_id int not null,
			gene_id	varchar(500),
			source	varchar(500),
			error_code int
);
create index idx_geneclusters_errorcode on gene_clusters(error_code);
create index idx_geneclusters_source on gene_clusters(source(1), error_code);
create index idx_geneclusters_clusterid on gene_clusters(gene_cluster_id);

drop table if exists cluster_summary;
create table if not exists cluster_summary(
		gene_cluster_id int not null,
		cap_gene_count int,
		cap_transcript_count int,
		vb_gene_count int,	
		vb_transcript_count int,
		cap_max_error int,
		vb_max_error int
);
create index idx_clustersum_clusterid on cluster_summary(gene_cluster_id);

drop table if exists transcript_mappings;
create table if not exists transcript_mappings(
			id       int not null AUTO_INCREMENT PRIMARY KEY,
			gene_cluster_id     int,
			cap_trans_id  varchar(500) default NULL,
			vb_trans_id   varchar(500) default NULL,
			map_type      varchar(500)
);
create index idx_trmap_cluster on transcript_mappings(gene_cluster_id);

drop table if exists gene_mappings;
create table if not exists gene_mappings(
			gene_cluster_id     int,
			cap_gene_id   varchar(500) default NULL,
			vb_gene_id    varchar(500) default NULL,
			map_type      varchar(500)
);
create index idx_gmappings_maptype on gene_mappings(map_type(6));
create index idx_gmappings_cap on gene_mappings(cap_gene_id, map_type(6));
create index idx_gmappings_vb on gene_mappings(vb_gene_id, map_type(6));

drop table if exists transcript_links;
create table if not exists transcript_links(
			id                  int not null AUTO_INCREMENT PRIMARY KEY,
			gene_cluster_id     int,
			cap_transcript_id   varchar(500),
			vb_transcript_id    varchar(500),
			link_group          varchar(500),
			link_rank           int,
			group_count		   int,
			link_status	       varchar(500)
						
);
create index idx_trlinks_cluster on transcript_links(gene_cluster_id, link_status(6));
create index idx_trlinks_cap on transcript_links(cap_transcript_id(6));
create index idx_trlinks_vb on transcript_links(vb_transcript_id(6));
create index idx_trlinks_status on transcript_links(link_status(3));

						
drop table if exists gene_events;
create table if not exists gene_events(
			id       int not null AUTO_INCREMENT PRIMARY KEY,
			vb_gene_id  TEXT default NULL,
			cap_gene_id TEXT default NULL,
			events  	varchar(500),
			vb_biotype varchar(500),
			cap_biotype varchar(500)
);


drop table if exists mapping_type;
create table if not exists mapping_type(
			map_type varchar(500),
			description varchar(200)
);
