#DROP TABLE snp;
#DROP TABLE reference;
#DROP TABLE population;
#DROP TABLE flags;

CREATE TABLE snp (

  rs                     VARCHAR(15) PRIMARY KEY,
  chr                    VARCHAR(6),
  pos			 INT NOT NULL,
  flags			 VARCHAR(500),
  ref_id		 INT NOT NULL,
  ref_base		 VARCHAR(1),
  alt_base		 VARCHAR(1),
  multi_genotypes        BOOLEAN,

  KEY pos_idx (ref_id,chr,pos)

);

CREATE TABLE reference (

  ref_id     INT NOT NULL AUTO_INCREMENT,
  name     	VARCHAR(10),
  alias         VARCHAR(500)
);

CREATE TABLE population (

  pop          VARCHAR(5) NOT NULL,
  rs           VARCHAR(15) NOT NULL,
  sample_size  INT,
  allele_freq  FLOAT,

  PRIMARY KEY (rs),
  KEY rs_pop_idx (rs, pop)
);

CREATE TABLE flags (
  id  VARCHAR(5) PRIMARY KEY,
  full VARCHAR(500)
);
