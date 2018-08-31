
-- ---------------------------------------------------------------------------------------------------

\set expected_version 94

\set ON_ERROR_STOP on

    -- warn that we detected the schema version mismatch:
SELECT ('The patch only applies to schema version '
    || CAST(:expected_version AS VARCHAR)
    || ', but the current schema version is '
    || meta_value
    || ', so skipping the rest.') as incompatible_msg
    FROM hive_meta WHERE meta_key='hive_sql_schema_version' AND meta_value!=CAST(:expected_version AS VARCHAR);

    -- cause division by zero only if current version differs from the expected one:
INSERT INTO hive_meta (meta_key, meta_value)
   SELECT 'this_should_never_be_inserted', 1 FROM hive_meta WHERE 1 != 1/CAST( (meta_key!='hive_sql_schema_version' OR meta_value=CAST(:expected_version AS VARCHAR)) AS INTEGER );

SELECT ('The patch seems to be compatible with schema version '
    || CAST(:expected_version AS VARCHAR)
    || ', applying the patch...') AS compatible_msg;


-- ----------------------------------<actual_patch> -------------------------------------------------


CREATE OR REPLACE VIEW resource_usage_stats AS
    SELECT a.logic_name || '(' || a.analysis_id || ')' analysis,
           w.meadow_type,
           rc.name || '(' || rc.resource_class_id || ')' resource_class,
           u.exit_status,
           count(*) workers,
           min(mem_megs) AS min_mem_megs, round(avg(mem_megs)*100)/100 AS avg_mem_megs, max(mem_megs) AS max_mem_megs,
           min(swap_megs) AS min_swap_megs, round(avg(swap_megs)*100)/100 AS avg_swap_megs, max(swap_megs) AS max_swap_megs,
           round(min(cpu_sec/lifespan_sec)*100)/100 AS min_cpu_perc, round(avg(cpu_sec/lifespan_sec)*100)/100 AS avg_cpu_perc, round(max(cpu_sec/lifespan_sec)*100)/100 AS max_cpu_perc
    FROM resource_class rc
    JOIN analysis_base a USING(resource_class_id)
    LEFT JOIN role r USING(analysis_id)
    LEFT JOIN worker w USING(worker_id)
    LEFT JOIN worker_resource_usage u USING (worker_id)
    GROUP BY a.analysis_id, w.meadow_type, rc.resource_class_id, u.exit_status
    ORDER BY a.analysis_id, w.meadow_type, rc.resource_class_id, u.exit_status;

CREATE OR REPLACE VIEW live_roles AS
    SELECT w.meadow_user, w.meadow_type, w.resource_class_id, rc.name resource_class_name, r.analysis_id, a.logic_name, count(*) AS num_workers
    FROM worker w
    JOIN role r USING(worker_id)
    LEFT JOIN resource_class rc ON w.resource_class_id=rc.resource_class_id
    LEFT JOIN analysis_base a USING(analysis_id)
    WHERE r.when_finished IS NULL
    GROUP BY w.meadow_user, w.meadow_type, w.resource_class_id, rc.name, r.analysis_id, a.logic_name;

CREATE OR REPLACE VIEW job_semaphores AS
    SELECT s.semaphore_id,
           fan_ab.logic_name AS controlling_logic_name,
           fan_j.job_id AS controlling_job_id,
           fan_j.status AS controlling_job_status,
           funnel_ab.logic_name AS semaphored_logic_name,
           s.dependent_job_id AS semaphored_job_id, s.dependent_semaphore_url
    FROM semaphore s
    LEFT JOIN job fan_j ON s.semaphore_id = fan_j.controlled_semaphore_id
    LEFT JOIN analysis_base fan_ab ON fan_ab.analysis_id = fan_j.analysis_id
    LEFT JOIN job funnel_j ON funnel_j.job_id = s.dependent_job_id
    LEFT JOIN analysis_base funnel_ab ON funnel_ab.analysis_id = funnel_j.analysis_id;

-- must drop then create to avoid errors due to typecasting
DROP VIEW IF EXISTS semaphore_summary;
CREATE VIEW semaphore_summary AS
   SELECT s.semaphore_id, COALESCE(funnel_ab.logic_name, '(see remote hive)') as semaphored_logic_name,
          s.dependent_job_id as semaphored_job_id, s.dependent_semaphore_url,
          COALESCE(fan_ab.logic_name, '(see remote hive)') as controlling_logic_name,
          COUNT(fan_j.job_id) as fan_jobs_this_analysis, s.remote_jobs_counter,
          COALESCE(to_char(min(fan_j.job_id), '999'), '(see remote hive)') as example_controlling_job_id
   FROM semaphore s
   LEFT JOIN job fan_j ON s.semaphore_id = fan_j.controlled_semaphore_id
   LEFT JOIN analysis_base fan_ab ON fan_ab.analysis_id = fan_j.analysis_id
   LEFT JOIN job funnel_j ON funnel_j.job_id = s.dependent_job_id
   LEFT JOIN analysis_base funnel_ab ON funnel_ab.analysis_id = funnel_j.analysis_id
   WHERE fan_j.status NOT IN ('DONE', 'PASSED_ON') OR fan_j.status IS NULL
   GROUP BY s.semaphore_id, funnel_ab.logic_name, s.dependent_job_id, fan_ab.logic_name, fan_j.analysis_id;


-- ----------------------------------</actual_patch> -------------------------------------------------


    -- increase the schema version by one and register the patch:
UPDATE hive_meta SET meta_value= (CAST(meta_value AS INTEGER) + 1) WHERE meta_key='hive_sql_schema_version';
INSERT INTO hive_meta (meta_key, meta_value) SELECT 'patched_to_' || meta_value, CURRENT_TIMESTAMP FROM hive_meta WHERE meta_key = 'hive_sql_schema_version';
