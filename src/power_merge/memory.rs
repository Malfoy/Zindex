//! Admission estimates, not a virtual-memory limit or an external sorter.
use anyhow::Result;
use std::{fs, path::Path, sync::mpsc, thread, time::Duration};

const GIB: u64 = 1 << 30;

pub(super) fn available_bytes() -> Option<u64> {
    let mut available = fs::read_to_string("/proc/meminfo").ok().and_then(|text| {
        text.lines().find_map(|line| {
            line.strip_prefix("MemAvailable:")?
                .split_whitespace()
                .next()?
                .parse::<u64>()
                .ok()
                .map(|n| n.saturating_mul(1024))
        })
    });
    if let (Some(n), Ok(arc)) = (
        available,
        fs::read_to_string("/proc/spl/kstat/zfs/arcstats"),
    ) {
        available = Some(n.saturating_add(reclaimable_arc(&arc)));
    }
    cap_available(available, cgroup_available_bytes())
}

fn cap_available(host: Option<u64>, cgroup: Option<u64>) -> Option<u64> {
    match (host, cgroup) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (a, b) => a.or(b),
    }
}

/// ARC is separate from the page cache already included in MemAvailable.
/// Only explicitly evictable ARC above its minimum is added, never all cache.
fn reclaimable_arc(text: &str) -> u64 {
    let value = |name: &str| {
        text.lines().find_map(|line| {
            let mut fields = line.split_whitespace();
            (fields.next()? == name)
                .then(|| fields.nth(1)?.parse::<u64>().ok())
                .flatten()
        })
    };
    let Some(size) = value("size") else {
        return 0;
    };
    let Some(minimum) = value("c_min") else {
        return 0;
    };
    let evictable = value("evictable").or_else(|| {
        Some(
            value("mru_evictable_data")?
                .saturating_add(value("mru_evictable_metadata")?)
                .saturating_add(value("mfu_evictable_data")?)
                .saturating_add(value("mfu_evictable_metadata")?),
        )
    });
    evictable.unwrap_or(0).min(size.saturating_sub(minimum))
}

fn cgroup_available_bytes() -> Option<u64> {
    let mut available: Option<u64> = None;
    // Respect visible cgroup-v2 ancestor limits as well as host RAM. Some
    // systems hide the cgroup filesystem; in that case retain MemAvailable.
    if let Ok(cgroup) = fs::read_to_string("/proc/self/cgroup") {
        if let Some(relative) = cgroup.lines().find_map(|s| s.strip_prefix("0::/")) {
            let root = Path::new("/sys/fs/cgroup");
            let mut path = root.join(relative);
            while path.starts_with(root) {
                let read = |name| {
                    fs::read_to_string(path.join(name))
                        .ok()?
                        .trim()
                        .parse::<u64>()
                        .ok()
                };
                if let (Some(limit), Some(current)) = (read("memory.max"), read("memory.current")) {
                    let free = limit.saturating_sub(current);
                    available = Some(available.map_or(free, |n| n.min(free)));
                }
                if path == root || !path.pop() {
                    break;
                }
            }
        }
    }
    available
}

pub(super) fn rss_bytes() -> Option<u64> {
    fs::read_to_string("/proc/self/status")
        .ok()?
        .lines()
        .find_map(|line| {
            line.strip_prefix("VmRSS:")?
                .split_whitespace()
                .next()?
                .parse::<u64>()
                .ok()
                .map(|n| n.saturating_mul(1024))
        })
}

pub(super) fn total_bytes(gb: Option<f64>) -> Option<u64> {
    gb.map(|n| (n * 1_000_000_000.0) as u64)
}

fn fixed_workspace(total: u64, baseline: u64, cgroup_free: Option<u64>) -> u64 {
    total
        .saturating_sub(baseline)
        .min(cgroup_free.unwrap_or(u64::MAX))
}

pub(super) fn workspace_budget(args: &super::MergeArgs) -> u64 {
    total_bytes(args.max_ram_gb)
        .map(|n| fixed_workspace(n, rss_bytes().unwrap_or(0), cgroup_available_bytes()).max(1))
        .unwrap_or_else(|| budget(args.build_memory_gib))
}

/// Active workspace estimates already include active resident arrays. In total
/// mode compare their reservations with total minus the phase's *initial*
/// retained baseline, never with free RAM that already excludes those arrays.
pub(super) fn total_admission_limit(total: u64, baseline: u64) -> u64 {
    fixed_workspace(
        total,
        baseline,
        cgroup_available_bytes()
            .map(|free| free.saturating_add(rss_bytes().unwrap_or(0).saturating_sub(baseline))),
    )
}

pub(super) struct RssGuard {
    stop: Option<mpsc::Sender<()>>,
    worker: Option<thread::JoinHandle<()>>,
    limit: Option<u64>,
}

impl RssGuard {
    pub(super) fn start(gb: Option<f64>) -> Result<Self> {
        let limit = total_bytes(gb);
        let mut guard = Self {
            stop: None,
            worker: None,
            limit,
        };
        guard.check()?;
        if let Some(limit) = limit {
            let (tx, rx) = mpsc::channel();
            guard.stop = Some(tx);
            eprintln!("total RAM policy: max_ram_bytes={limit}; decimal GB; fixed total RSS; cgroups respected; poll_interval_ms=20; transient overshoot between polls is possible");
            guard.worker = Some(thread::spawn(move || {
                while rx.recv_timeout(Duration::from_millis(20)).is_err() {
                    if let Err(error) = check_rss(limit) {
                        eprintln!("{error:#}");
                        // Abort all Rayon work immediately. No partial index is
                        // reported as complete; staging is left for inspection.
                        std::process::exit(86);
                    }
                }
            }));
        }
        Ok(guard)
    }
    pub(super) fn check(&self) -> Result<()> {
        if let Some(limit) = self.limit {
            check_rss(limit)?;
        }
        Ok(())
    }
}

fn check_rss(limit: u64) -> Result<()> {
    let rss =
        rss_bytes().ok_or_else(|| anyhow::anyhow!("max RAM guard cannot read process RSS"))?;
    if rss > limit {
        eprintln!(
            "memory-limit-exceeded\t{}",
            serde_json::json!({
            "kind":"actual_rss_limit", "rss_bytes":rss, "limit_bytes":limit,
            "poll_interval_ms":20, "pid":std::process::id()})
        );
    }
    anyhow::ensure!(rss <= limit, "max RAM exceeded: actual process RSS {rss} bytes > --max-ram-gb ceiling {limit} bytes (measured RSS, not address space)");
    if let Some(free) = cgroup_available_bytes() {
        if free == 0 {
            eprintln!(
                "memory-limit-exceeded\t{}",
                serde_json::json!({
                "kind":"actual_cgroup_limit", "rss_bytes":rss, "limit_bytes":limit,
                "cgroup_available_bytes":free,"poll_interval_ms":20,"pid":std::process::id()})
            );
        }
        anyhow::ensure!(free > 0, "max RAM guard: actual cgroup memory capacity exhausted; RSS={rss} bytes, requested ceiling={limit} bytes");
    }
    Ok(())
}

impl Drop for RssGuard {
    fn drop(&mut self) {
        if let Some(stop) = self.stop.take() {
            let _ = stop.send(());
        }
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

pub(super) fn budget(gib: Option<f64>) -> u64 {
    gib.map(|n| (n * GIB as f64) as u64)
        .unwrap_or_else(|| available_bytes().unwrap_or(8 * GIB).saturating_mul(4) / 5)
        .max(1)
}

fn limit_for_available(budget: u64, explicit: bool, available: Option<u64>) -> u64 {
    available
        .map(|n| budget.min(if explicit { n } else { n.saturating_mul(4) / 5 }))
        .unwrap_or(budget)
}

/// An explicit workspace budget already expresses the caller's headroom choice.
/// Do not silently reduce it to the automatic80% policy a second time.
pub(super) fn admission_limit(budget: u64, explicit: bool) -> u64 {
    limit_for_available(budget, explicit, available_bytes())
}

pub(super) fn input_estimate(path: &Path) -> Result<u64> {
    // Extraction stores occurrences before deduplication. Compressed size is
    // only a heuristic (especially for very repetitive data), never a bound.
    let bytes = fs::metadata(path)?.len();
    Ok(bytes
        .saturating_mul(if crate::is_compressed_path(path) {
            64
        } else {
            24
        })
        .max(8 << 20))
}

/// Pages written by the u64 occurrence vector before duplicate removal.
/// Unused Vec capacity reserves address space, not resident RAM.
pub(super) fn feature_resident_bytes(occurrences: u64) -> u64 {
    occurrences.saturating_mul(8)
}

/// Opt-in construction diagnostics. No probes or timers run when disabled.
pub fn trace_memory(stage: &str, detail: serde_json::Value) {
    if std::env::var_os("ZORINDEX_MEMORY_TRACE").is_none() {
        return;
    }
    let status = fs::read_to_string("/proc/self/status").unwrap_or_default();
    let read = |name: &str| {
        status
            .lines()
            .find_map(|s| {
                s.strip_prefix(name)?
                    .split_whitespace()
                    .next()?
                    .parse::<u64>()
                    .ok()
            })
            .map(|n| n.saturating_mul(1024))
    };
    eprintln!(
        "memory-trace\t{}",
        serde_json::json!({"stage": stage,
        "rss_bytes": read("VmRSS:"), "hwm_bytes": read("VmHWM:"), "detail": detail})
    );
}

fn compact_graph(features: usize, slots: usize, hashes: usize) -> bool {
    features <= u32::MAX as usize
        && slots <= u32::MAX as usize
        && features
            .checked_mul(hashes)
            .is_some_and(|n| n <= u32::MAX as usize)
}

pub(super) fn prediction_kind(features: usize, slots: usize, hashes: usize) -> &'static str {
    if compact_graph(features, slots, hashes) {
        "compact-resident-v2"
    } else {
        "wide-capacity-bound-unvalidated"
    }
}

/// Expected graph/assignment RSS, without the admission safety margin.
/// Compact-graph residency is validated against real multi-Gb constructions.
/// Wide graphs keep their explicit conservative fallback until measured.
pub(super) fn group_prediction(
    features: usize,
    slots: usize,
    hashes: usize,
    width: usize,
    input_peak: u64,
) -> u64 {
    if !compact_graph(features, slots, hashes) {
        return capacity_bound(features, slots, hashes, width, input_peak);
    }
    let n = features as u64;
    let s = slots as u64;
    let keys = input_peak.max(n.saturating_mul(8));
    let load = (features as f64) * (hashes as f64) / (slots.max(1) as f64);
    let empty = (-load).exp();
    // The u32 location array is fully initialized. Heap nodes are written only
    // for degree>=2 cells: under uniform hashing, P(degree>=2)=1-e^-λ(1+λ).
    // The heap never acquires new cells as degrees only decrease.
    let occupied = (1.0 - empty * (1.0 + load)).clamp(0.0, 1.0);
    let heap = (s as f64) * (4.0 + 8.0 * occupied);
    // Geometric VecDeque capacity for initially degree-one cells. This is a
    // prediction, not a bound on every subsequent peeling wave.
    let initial_one = ((s as f64) * load * empty).ceil() as u64;
    let queue = if initial_one == 0 {
        0
    } else {
        initial_one
            .checked_next_power_of_two()
            .unwrap_or(u64::MAX)
            .saturating_mul(4)
    };
    // Keys, payloads, both incidence arrays, active flags and peel-order data.
    let graph = keys
        .saturating_add(n.saturating_mul(8 + 8 * hashes as u64 + 1 + 8))
        .saturating_add(s.saturating_mul(12))
        .saturating_add(heap.ceil() as u64)
        .saturating_add(queue);
    // Assignment runs after adjacency, degrees, queue and heap have been freed.
    let assignment = keys
        .saturating_add(n.saturating_mul(8 + 4 * hashes as u64 + 8))
        .saturating_add(
            s.saturating_mul((width.div_ceil(8).next_power_of_two() + width.div_ceil(8)) as u64),
        );
    graph.max(assignment).max(8 << 20)
}

/// Admission uses the RSS prediction plus a separately stated 10% margin.
pub(super) fn group_estimate(
    features: usize,
    slots: usize,
    hashes: usize,
    width: usize,
    input_peak: u64,
) -> u64 {
    let peak = group_prediction(features, slots, hashes, width, input_peak);
    peak.saturating_add(peak.div_ceil(10))
}

/// Conservative capacity-based fallback; deliberately not an expected RSS.
pub(super) fn capacity_bound(
    features: usize,
    slots: usize,
    hashes: usize,
    width: usize,
    input_peak: u64,
) -> u64 {
    let narrow = features <= u32::MAX as usize
        && slots <= u32::MAX as usize
        && features
            .checked_mul(hashes)
            .is_some_and(|n| n <= u32::MAX as usize);
    let id = if narrow { 4u64 } else { 8 };
    let n = features as u64;
    let s = slots as u64;
    // Retained occurrence-vector capacity survives deduplication. Merged keys
    // and payloads may have up to 2x spare capacity; never charge file size.
    let keys = input_peak
        .checked_next_power_of_two()
        .unwrap_or(u64::MAX)
        .max(n.saturating_mul(16));
    // The compact indexed heap uses 8-byte nodes and 4-byte locations per cell.
    // For the wide lazy heap, each initial/update push accounts for an incidence:
    // initial(degree>1) + updates(new_degree>1) = max(degree-1,0).
    // In no-heuristic mode, capacity cannot exceed pow2(max(S,h*N)). Other
    // heuristics only reinsert popped entries, so do not increase this bound.
    let heap = if narrow {
        s.saturating_mul(12)
    } else {
        s.max(n.saturating_mul(hashes as u64))
            .checked_next_power_of_two()
            .unwrap_or(u64::MAX)
            .saturating_mul(16)
    };
    let graph = keys
        .saturating_add(n.saturating_mul(16 + 2 * hashes as u64 * id + 1 + 2 * id))
        .saturating_add(s.saturating_mul(4 + 4 * id))
        .saturating_add(heap);
    // Adjacency/degrees/heap have been freed before cell assignment. Bound
    // that separate phase too rather than summing non-overlapping phases.
    let assignment = keys
        .saturating_add(n.saturating_mul(16 + hashes as u64 * id + 2 * id))
        .saturating_add(
            s.saturating_mul((width.div_ceil(8).next_power_of_two() + width.div_ceil(8)) as u64),
        );
    let peak = graph.max(assignment);
    peak.max(8 << 20)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(target_os = "linux")]
    #[test]
    fn watchdog_stops_measured_growth_in_child_process() {
        const CHILD: &str = "ZINDEX_TEST_RSS_WATCHDOG_CHILD";
        if std::env::var_os(CHILD).is_some() {
            let baseline = rss_bytes().expect("Linux process RSS");
            let limit = baseline + (8 << 20);
            let _guard = RssGuard::start(Some(limit as f64 / 1_000_000_000.0)).unwrap();
            // Touch and retain resident pages after starting the real watchdog.
            // No construction admission estimate participates in this fixture.
            let mut pages = vec![0x5au8; 16 << 20];
            std::hint::black_box(&mut pages);
            thread::sleep(Duration::from_secs(5));
            std::hint::black_box(&pages);
            panic!("watchdog did not stop a process exceeding its RSS ceiling");
        }

        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "power_merge::memory::tests::watchdog_stops_measured_growth_in_child_process",
                "--nocapture",
            ])
            .env(CHILD, "1")
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(86), "{stderr}");
        let evidence = stderr
            .lines()
            .find_map(|line| line.strip_prefix("memory-limit-exceeded\t"))
            .expect("watchdog must report measured RSS evidence");
        let evidence: serde_json::Value = serde_json::from_str(evidence).unwrap();
        assert_eq!(evidence["kind"], "actual_rss_limit");
        assert!(
            evidence["rss_bytes"].as_u64().unwrap() > evidence["limit_bytes"].as_u64().unwrap()
        );
    }

    #[test]
    fn total_decimal_budget_has_no_memavailable_clamp_or_double_counting() {
        assert_eq!(total_bytes(Some(500.0)), Some(500_000_000_000));
        assert_eq!(
            fixed_workspace(500_000_000_000, 20_000_000_000, None),
            480_000_000_000
        );
        // Host MemAvailable=206GB deliberately has no place in fixed mode.
        assert_eq!(
            fixed_workspace(500_000_000_000, 20_000_000_000, Some(170_000_000_000)),
            170_000_000_000
        );
        assert_eq!(fixed_workspace(10, 20, None), 0);
    }
    #[test]
    fn arc_reclaimability_is_explicit_bounded_and_then_cgroup_capped() {
        let arc = "size 4 205000000000\nc_min 4 12000000000\nmru_evictable_data 4 120000000000\nmru_evictable_metadata 4 1000000000\nmfu_evictable_data 4 71000000000\nmfu_evictable_metadata 4 1000000000\n";
        assert_eq!(reclaimable_arc(arc), 193_000_000_000);
        assert_eq!(
            cap_available(
                Some(206_000_000_000 + reclaimable_arc(arc)),
                Some(250_000_000_000)
            ),
            Some(250_000_000_000)
        );
        assert_eq!(
            reclaimable_arc("size 4 100\nc_min 4 20\nevictable 4 500\n"),
            80
        );
        assert_eq!(
            reclaimable_arc("size 4 100\nc_min 4 120\nevictable 4 500\n"),
            0
        );
        assert_eq!(reclaimable_arc("size 4 100\nc_min 4 20\n"), 0);
    }
    #[test]
    fn explicit_budget_is_not_reduced_by_automatic_headroom() {
        let gb = 1_000_000_000;
        assert_eq!(
            limit_for_available(500 * gb, true, Some(510 * gb)),
            500 * gb
        );
        assert_eq!(
            limit_for_available(500 * gb, false, Some(510 * gb)),
            408 * gb
        );
        assert_eq!(
            limit_for_available(500 * gb, true, Some(250 * gb)),
            250 * gb
        );
        assert_eq!(limit_for_available(500 * gb, true, None), 500 * gb);
    }
    #[test]
    fn estimates_account_for_wide_graph_and_do_not_overflow() {
        let narrow = group_estimate(1000, 2048, 8, 8, 0);
        let wide = group_estimate(536_870_912, 1 << 30, 8, 8, 0);
        assert!(wide > narrow);
        assert!(wide > 536_870_912u64 * 128);
        assert_eq!(
            group_estimate(usize::MAX, usize::MAX, 32, 64, u64::MAX),
            u64::MAX
        );
        assert_eq!(budget(Some(0.5)), GIB / 2);
    }

    #[test]
    fn actual_feature_reservation_and_hash_width_control_graph_estimate() {
        assert_eq!(feature_resident_bytes(3_000_000_000), 24_000_000_000);
        let n = 1_000_000_000;
        let h4 = group_estimate(n, 1 << 30, 4, 8, 32 * GIB);
        let h8 = group_estimate(n, 1 << 30, 8, 8, 32 * GIB);
        assert!(h4 < 150 * GIB);
        assert!(h8 > 2 * h4);
        assert!(h8 < 400 * GIB);
    }

    #[test]
    fn resident_prediction_matches_completed_large_genome_constructions() {
        // BWT native OCM31/window49, fp8, two-pass singleton constructions.
        // r10-h4 was derivation; h6/h8 predictions were registered before those
        // builds completed. r09 is a separate genome. No fitted constants here.
        let observed: [(usize, u64, usize, u64); 5] = [
            (203_688_336, 785_004_796, 4, 22_562_078_720u64),
            (203_688_336, 785_004_796, 6, 25_910_657_024),
            (203_688_336, 785_004_796, 8, 29_209_948_160),
            (274_986_295, 758_953_387, 4, 31_824_310_272),
            (274_986_295, 758_953_387, 6, 36_565_446_656),
        ];
        for (n, o, h, rss) in observed {
            let slots = n.next_power_of_two();
            let p = group_prediction(n, slots, h, 8, feature_resident_bytes(o));
            let error = (p as f64 / rss as f64 - 1.0).abs();
            assert!(error < 0.01, "h={h},N={n}: predicted={p},measured={rss}");
            assert!(group_estimate(n, slots, h, 8, feature_resident_bytes(o)) > rss);
        }
        assert_eq!(
            prediction_kind(1_000_000_000, 1 << 30, 4),
            "compact-resident-v2"
        );
        assert_eq!(
            prediction_kind(1_000_000_000, 1 << 30, 8),
            "wide-capacity-bound-unvalidated"
        );
    }
}
