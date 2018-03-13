drop table if exists exon;
create table if not exists exon(
			exon_id  varchar(500),
			source   varchar(500),
			scaffold varchar(500),
			strand   varchar(500),
			start    int,
			end      int			
);

drop table if exists cds;
create table if not exists cds(
			cds_parent_id varchar(500),
			start         int,
			end           int

);

drop table if exists gene_model;
create table if not exists gene_model(
			id 			  int not null AUTO_INCREMENT PRIMARY KEY,
			exon_id       varchar(500),
			transcript_id varchar(500),
			gene_id       varchar(500),
			source        varchar(500)
);


drop table if exists exon_mappings;
create table if not exists exon_mappings(
			id       int not null AUTO_INCREMENT PRIMARY KEY,
			cap_exon_id   varchar(500),
			vb_exon_id    varchar(500),
			map_type      varchar(500)
			
);

drop table if exists gene_clusters;
create table if not exists gene_clusters(
			gene_cluster_id int not null,
			gene_id	varchar(500),
			source	varchar(500)
);

drop table if exists cluster_summary;
create table if not exists cluster_summary(
		gene_cluster_id int not null,
		cap_gene_count int,
		cap_transcript_count int,
		vb_gene_count int,	
		vb_transcript_count int
);

drop table if exists transcript_mappings;
create table if not exists transcript_mappings(
			id       int not null AUTO_INCREMENT PRIMARY KEY,
			gene_cluster_id     int,
			cap_trans_id  varchar(500) default NULL,
			vb_trans_id   varchar(500) default NULL,
			map_type      varchar(500)
);

drop table if exists gene_mappings;
create table if not exists gene_mappings(
			gene_cluster_id     int,
			cap_gene_id   varchar(500) default NULL,
			vb_gene_id    varchar(500) default NULL,
			map_type      varchar(500)
);


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

						
drop table if exists gene_events;
create table if not exists gene_events(
			id       int not null AUTO_INCREMENT PRIMARY KEY,
			vb_gene_id  TEXT default NULL,
			cap_gene_id TEXT default NULL,
			events  	varchar(500)   
);


drop table if exists mapping_type;
create table if not exists mapping_type(
			map_type varchar(500),
			description varchar(200)
);
