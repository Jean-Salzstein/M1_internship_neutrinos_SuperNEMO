Ce dossier GITHUB a été créé à l'occasion de mon stage de Master 1, pour rendre public mon algorithme de sélection des traces à partir notamment d'une clusterisation basique de SuperNEMO. Il contient l'événement 2589 de la série 1307 de mesures SuperNEMO (3h, source neutrons AmBe).

Exemple d'utilisation :

```
root [0] .L Code_GITHUB.cpp
root [1] show_clusterised_event()
=> run file "Classification_events_1307_github_version.root" opened (1 events)
=> run file "snemo_run-1307_event-2589_github_version.root" opened (1 events)
Number of : e - 2 - g - 1 - alpha - 0 - weirdies - 0.
```
![Top_view_2589_internal](https://github.com/Jean-Salzstein/M1_internship_neutrinos_SuperNEMO/assets/173078447/ec4b1b6b-06b0-42b1-aa2b-ce4b282ffcd6)
![Calo_walls_view_2589_internal](https://github.com/Jean-Salzstein/M1_internship_neutrinos_SuperNEMO/assets/173078447/29901661-4758-4221-9227-0748f34d590a)

```
root [2] cout << "Nombre d'electrons : " << nb_total_elec << endl << "Correlation temporelle entre les modules optiques : " << dt_good_corr->at(0) << "ns." << endl;
Nombre d'electrons : 2
Correlation temporelle entre les modules optiques : 1.83413ns.
```

