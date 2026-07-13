#include "NuHepMC/HepMC3Features.hxx"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "NuHepMC/EventUtils.hxx"
#include "NuHepMC/FATXUtils.hxx"
#include "NuHepMC/ReaderUtils.hxx"
#include "NuHepMC/Reader.hxx"

#include <iostream>
#include <sstream>
#include <vector>

using namespace NuHepMC;

struct ParticleRecord {
  int id;       // HepMC internal particle id
  int pid;      // PDG code
  int status;   // HepMC / NuHepMC status
  double px;
  double py;
  double pz;
  double E;

  int production_vertex_id = 0; // 0 or -1 sentinel if none
  int end_vertex_id = 0;        // 0 or -1 sentinel if none

  bool is_beam = false;
  bool is_target = false;
  bool is_primary_process = false;
  bool is_final_physical = false;
  bool is_remnant_pseudoparticle = false;
  double p() const {
    return std::sqrt(px*px + py*py + pz*pz);
  }
};

struct EventRecord {
  int event_number;
  int process_id;
  std::vector<ParticleRecord> particles;
};

ParticleRecord MakeParticleRecord(HepMC3::ConstGenParticlePtr p) {
  auto mom = p->momentum();
  ParticleRecord rec;
  rec.id = p->id();
  rec.pid = p->pid();
  rec.status = p->status();
  rec.px = mom.px();
  rec.py = mom.py();
  rec.pz = mom.pz();
  rec.E  = mom.e();

  rec.production_vertex_id =
      p->production_vertex() ? p->production_vertex()->id() : -1;
  rec.end_vertex_id =
      p->end_vertex() ? p->end_vertex()->id() : -1;

  rec.is_beam = (p->status() == 4);
  rec.is_target = (p->status() == 20);
  rec.is_primary_process = (p->status() == 3);
  rec.is_final_physical = (p->status() == 1);
  rec.is_remnant_pseudoparticle = (p->pid() == 2009900000);

  return rec;
}

EventRecord MakeEventRecord(HepMC3::GenEvent &evt) {
  EventRecord rec;
  rec.event_number = evt.event_number();
  rec.process_id = ER3::ReadProcessID(evt);

  rec.particles.reserve(evt.particles().size());
  for (auto const &p : evt.particles()) {
    rec.particles.push_back(MakeParticleRecord(p));
  }
  return rec;
}

void PrintParticle(const ParticleRecord &p, size_t idx = 0) {
  if(p.is_primary_process || p.is_remnant_pseudoparticle){
    return;
  }
  std::cout << "    [" << idx << "] "
            << "id=" << p.id
            << ", pid=" << p.pid
            << ", status=" << p.status
            << ", p=("
            << p.px << ", "
            << p.py << ", "
            << p.pz << ", E=" << p.E
            << ") MeV"
            << std::endl;
}

void PrintEvent(const EventRecord &ev) {
  std::cout << "====================================" << std::endl;
  std::cout << "Event: " << ev.event_number << std::endl;
  std::cout << "Process ID: " << ev.process_id << std::endl;

  std::cout << "N particles: " << ev.particles.size() << std::endl;
  for (size_t i = 0; i < ev.particles.size(); ++i) {
    PrintParticle(ev.particles[i], i);
  }
  std::cout << "====================================" << std::endl;
}

void GetMomentum(HepMC3::GenEvent &evt){
  std::cout << "Evt: " << evt.event_number() << std::endl;
  for (auto p : evt.particles()) {
    auto mom = p->momentum();
    std::cout << "PID: " << p->pid() << std::endl;
    std::cout << mom.px() << " "
              << mom.py() << " "
              << mom.pz() << " "
              << mom.e()  << std::endl;
    if(p->pid() == 13){
      double px = mom.px();
      double py = mom.py();
      double pz = mom.pz();
      double p_mag = sqrt(px*px + py*py + pz*pz);
      std::cout<< p_mag << std::endl;
    }
  }
}

int main(int argc, char const *argv[]) {
  std::vector<EventRecord> events;

  if (argc < 3) {
    std::cout << "[RUNLIKE]: " << argv[0] << " <infile.hepmc3> <event index>" << std::endl;
    return 1;
  }

  std::string inf = argv[1];
  int event_id = std::stoi(argv[2]);

  auto rdr = std::make_unique<NuHepMC::Reader>(inf);
  if (!rdr) {
    std::cout << "Failed to instantiate HepMC3::Reader from " << inf << std::endl;
    return 1;
  }
  HepMC3::GenEvent evt;
  rdr->read_event(evt);
  if (rdr->failed()) {
    std::cout << "Failed to read the first event from " << inf << std::endl;
    return 1;
  }
  // grab a copy of the GenRunInfo
  auto run_info = rdr->run_info();
  if (!run_info) {
    std::cout << "Warning, after reading the first event, we have no "
                 "GenRunInfo."
              << std::endl;
  }

  // re-open the file to start from the first event
  rdr = std::make_unique<NuHepMC::Reader>(inf);
  while (true) {

    // read an event and check that you haven't finished the file.
    // A file only knows it has read the last event after trying to read the
    // next one
    rdr->read_event(evt);
    if (rdr->failed()) {
      break;
    }
    // process event: put whatever you like here!
    // GetMomentum(evt);
    events.push_back(MakeEventRecord(evt));
  }

  PrintEvent(events[event_id]);
}
