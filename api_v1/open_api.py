import json
from datetime import datetime
import dateutil
from flask import Blueprint, request, Response
from schema.models import Enum, date, Ontology, ErrorResponse, SpeciesResponse, Dataset, DatasetsResponse, SubjectDataset, SubjectDatasetResponse, Sample, \
    SampleMetadataResponse, Repertoire, DataProcessing, SampleProcessing, CellProcessing, NucleicAcidProcessing, SequencingRun, LibraryGenerationMethod, TemplateClass, \
    CompleteSequences, PhysicalLinkage, SequencingData, FileTypeEnum, ReadDirectionEnum, PairedReadDirectionEnum, Subject, SexEnum, SubjectGenotype, GenotypeSet, Diagnosis, \
    Study, KeywordsStudyEnum, GenotypeClassListItem, AllSubjectsGenotypeResponse, AllSamplesMetadataResponse, PCRTarget
from pydantic.fields import FieldInfo
from pydantic import BaseModel
from typing import Any, Union, get_args, get_origin
from api.system.system import digby_protected
from api.genomic import genomic
from api.vdjbase import vdjbase
from flask_restx import Resource
from app import vdjbase_dbs, genomic_dbs


api_bp = Blueprint('api_v1', __name__)


def custom_jsonify(obj):
    """Custom JSON encoder for special object types."""

    def encode_obj(o):
        if isinstance(o, Enum):
            return o.value
        if isinstance(o, BaseModel):
            return o.model_dump()
        if isinstance(o, date):
            return o.isoformat()
        if isinstance(o, dict):
            return {k: encode_obj(v) for k, v in o.items()}
        if isinstance(o, list):
            return [encode_obj(i) for i in o]
        return o
    
    return Response(
        json.dumps(encode_obj(obj), indent=4),
        mimetype='application/json'
    )


def common_lookup(binomial):
    for lookup_db in [genomic_dbs, vdjbase_dbs]:
        for species, datasets in lookup_db.items():
            for ds_name, ds_data in datasets.items():
                if ds_data.binomial == binomial:
                    return species
    return None


@api_bp.route('/<type>/species', methods=['GET'])
@digby_protected()
def get_species(type):
    """Get species list based on type."""
    if type not in ['genomic', 'airrseq']:
        error_response = ErrorResponse(message="dataset type not valid")
        return error_response.model_dump_json(), 500

    species_list = []
    ontology_list = []

    lookup_dbs = genomic_dbs if type == "genomic" else vdjbase_dbs
    for sp, datasets in lookup_dbs.items():
        for ds_name, ds_data in datasets.items():
            if sp not in species_list:
                species_list.append(sp)
                ontology_obj = Ontology(id=ds_data.taxid, label=ds_data.binomial)
                ontology_list.append(ontology_obj)

    species_response_obj = SpeciesResponse(species=ontology_list)

    try:
        return species_response_obj.model_dump_json(), 200
    except Exception as e:
        error_response = ErrorResponse(message=str(e))
        return error_response.model_dump_json(), 500


@api_bp.route('/<type>/datasets/<species>', methods=['GET'])
@digby_protected()
def get_species_datasets(type, species):
    """Get datasets for a species based on type."""
    if type not in ['genomic', 'airrseq']:
        error_response = ErrorResponse(message="dataset type not valid")
        return error_response.model_dump_json(), 500

    dataset_list = []
    lookup_dbs = genomic_dbs if type == "genomic" else vdjbase_dbs

    for _, datasets in lookup_dbs.items():
        for ds_name, ds_data in datasets.items():
            if ds_data.binomial == species:
                locus = ds_name
                dataset_obj = Dataset(dataset=locus, locus=locus, type=type, revision_date=ds_data.created)
                dataset_list.append(dataset_obj)

    dataset_response = DatasetsResponse(datasets=dataset_list)
    try:
        return dataset_response.model_dump_json(), 200

    except Exception as e:
        error_response = ErrorResponse(message=str(e))
        return error_response.model_dump_json(), 500


@api_bp.route('/<type>/subjects/<species>/<dataset>', methods=['GET'])
@digby_protected()
def get_subject_datasets(type, species, dataset):
    species = common_lookup(species)

    if not species:
        error_response = ErrorResponse(message="species not found")
        return error_response.model_dump_json(), 500

    """Get subject datasets for a species and dataset based on type."""
    if type == "genomic":
        try:
            subjects_list = genomic.SubjectsAPI(Resource)
            subjects_list = subjects_list.get(species, dataset)[0]
            dataset_list = []
            for sample in subjects_list.get('samples'):
                subject_identifier = sample['sample_name'].split('_')[0:1]
                subject_identifier = '_'.join(sample['sample_name'].rsplit('_', 1)[:-1])

                subject_dataset_obj = SubjectDataset(id=sample['sample_id'],
                                                     study_name=sample['study_name'],
                                                     subject_identifier=subject_identifier,
                                                     sample_identifier=sample['sample_name'],
                                                     dataset=sample['dataset'])
                dataset_list.append(subject_dataset_obj)
            subject_dataset_response_obj = SubjectDatasetResponse(subject_datasets=dataset_list)
          
            try:
                return subject_dataset_response_obj.model_dump_json(), 200

            except Exception as e:
                error_response = ErrorResponse(message=str(e))
                return error_response.model_dump_json(), 500

        except Exception as e:
            error_response = ErrorResponse(message=str(e))
            return error_response.model_dump_json(), 500

    elif type == "airrseq":
        try:
            subjects_list = vdjbase.SamplesApi(Resource)
            request.args = {'cols': '["sample_name", "study_id", "subject_id"]'}
            subjects_list = subjects_list.get(species, dataset)[0]
            dataset_list = []
            for sample in subjects_list.get('samples'):
                subject_dataset_obj = SubjectDataset(id=sample['sample_name'], 
                                                     study_name=sample['study_id'], 
                                                     subject_identifier=sample['subject_id'], 
                                                     sample_identifier=sample['sample_id'],
                                                     dataset=sample['dataset'])
                dataset_list.append(subject_dataset_obj)
            subject_dataset_response_obj = SubjectDatasetResponse(subject_datasets=dataset_list)
          
            try:
                return subject_dataset_response_obj.model_dump_json(), 200

            except Exception as e:
                error_response = ErrorResponse(message=str(e))
                return error_response.model_dump_json(), 500
        
        except Exception as e:
            error_response = ErrorResponse(message=str(e))
            return error_response.model_dump_json(), 500


@api_bp.route('/<type>/subject_genotype/<species>/<subject>', methods=['GET'])
@digby_protected()
def get_sample_genotype(type, species, subject):
    species = common_lookup(species)

    if not species:
        error_response = ErrorResponse(message="species not found")
        return error_response.model_dump_json(), 400

    if type not in ['genomic', 'airrseq']:
        error_response = ErrorResponse(message="dataset type not valid")
        return error_response.model_dump_json(), 400

    if type == "genomic":
        genotype_list = genomic.GenotypeApi(Resource)
        genotype_list = genotype_list.get(species, subject)

    elif type == "airrseq":
        genotype_list = vdjbase.GenotypeApi(Resource)
        genotype_list = genotype_list.get(species, subject)
    
    if not genotype_list or not genotype_list[0]:
        error_response = ErrorResponse(message="Subject not found")
        return error_response.model_dump_json(), 400
    
    set = genotype_list[0]['GenotypeSet']

    genotype_set = GenotypeSet(**set)
    return genotype_set.model_dump_json(), 200


@api_bp.route('/<type>/all_subjects_genotype/<species>', methods=['GET'])
@digby_protected()
def get_all_subjects_genotype(type, species):
    species = common_lookup(species)

    if not species:
        error_response = ErrorResponse(message="species not found")
        return error_response.model_dump_json(), 400

    if type not in ['genomic', 'airrseq']:
        error_response = ErrorResponse(message="dataset type not valid")
        return error_response.model_dump_json(), 400

    if type == "genomic":
        genotype_list = genomic.AllSubjectsGenotypeApi(Resource)
        genotype_list = genotype_list.get(species)

    elif type == "airrseq":
        genotype_list = vdjbase.AllSubjectsGenotypeApi(Resource)
        genotype_list = genotype_list.get(species)
    
    if not genotype_list or not genotype_list[0]:
        error_response = ErrorResponse(message="Subject not found")
        return error_response.model_dump_json(), 400
    
    sets = genotype_list[0]
    sets = [GenotypeClassListItem(subject_name=s['subject_name'], genotypeSet=s['GenotypeSet']) for s in sets]

    return AllSubjectsGenotypeResponse(genotype_class_list=sets).model_dump_json(), 200


@api_bp.route('/<type>/sample_metadata/<species>/<dataset>/<sample>', methods=['GET'])
@digby_protected()
def get_sample_metadata(type, species, dataset, sample):
    """Get metadata for a specific sample."""
    species = common_lookup(species)

    if not species:
        error_response = ErrorResponse(message="species not found")
        return error_response.model_dump_json(), 400

    if type == "genomic":
        try:
            sample_info = genomic.SampleInfoApi(Resource)
            sample_info = sample_info.get(species, dataset, sample)
            if not sample_info or not sample_info[0]:
                error_response = ErrorResponse(message="Sample not found")
                return error_response.model_dump_json(), 400

            rep_obj = SampleMetadataResponse(Repertoire=create_repertoire_obj(sample_info[0]))
            return custom_jsonify(rep_obj.model_dump(by_alias=True)), 200

        except Exception as e:
            error_response = ErrorResponse(message=str(e))
            return error_response.model_dump_json(), 500

    elif type == "airrseq":
        try:
            sample_info = vdjbase.SampleInfoApi(Resource)
            sample_info = sample_info.get(species, dataset, sample)
            if not sample_info or not sample_info[0]:
                error_response = ErrorResponse(message="Sample not found")
                return error_response.model_dump_json(), 400

            rep_obj = SampleMetadataResponse(Repertoire=create_repertoire_obj(sample_info[0]))
            return custom_jsonify(rep_obj.model_dump(by_alias=True)), 200

        except Exception as e:
            error_response = ErrorResponse(message=str(e))
            return error_response.model_dump_json(), 500
    else:
        error_response = ErrorResponse(message=str("No such type"))
        return error_response.model_dump_json(), 500


@api_bp.route('/<type>/all_samples_metadata/<species>/<dataset>', methods=['GET'])
@digby_protected()
def get_all_samples_metadata(type, species, dataset):
    """Get metadata all samples in a dataset."""
    species = common_lookup(species)

    if not species:
        error_response = ErrorResponse(message="species not found")
        return error_response.model_dump_json(), 400

    if type == "genomic":
        try:
            sample_info = genomic.AllSamplesInfoApi(Resource)
            sample_info = sample_info.get(species, dataset)
            if not sample_info or not sample_info[0]:
                error_response = ErrorResponse(message="Sample not found")
                return error_response.model_dump_json(), 400

        except Exception as e:
            error_response = ErrorResponse(message=str(e))
            return error_response.model_dump_json(), 500

    elif type == "airrseq":
        try:
            sample_info = vdjbase.AllSamplesInfoApi(Resource)
            sample_info = sample_info.get(species, dataset)
            if not sample_info or not sample_info[0]:
                error_response = ErrorResponse(message="Sample not found")
                return error_response.model_dump_json(), 400

        except Exception as e:
            error_response = ErrorResponse(message=str(e))
            return error_response.model_dump_json(), 500
    else:
        error_response = ErrorResponse(message=str("No such type"))
        return error_response.model_dump_json(), 500

    repertoires = [create_repertoire_obj(s) for s in sample_info[0]]
    rep_obj = AllSamplesMetadataResponse(Repertoire=repertoires)
    return custom_jsonify(rep_obj.model_dump(by_alias=True)), 200


def create_repertoire_obj(subject_info):
    """Create a Repertoire object from subject information."""
    subject_info = fill_missing_required_fields(Repertoire,  subject_info)

    repertoire_id = f"{subject_info['patient_name']}_{subject_info['study_id']}_{subject_info['pcr_target_locus']}"
    if subject_info['repertoire_id']:
        repertoire_id += '_' + subject_info['repertoire_id']

    rep_object = Repertoire(repertoire_id=repertoire_id, 
                            repertoire_name=subject_info["repertoire_name"],
                            repertoire_description=subject_info["repertoire_description"],
                            study=create_study_object(subject_info),
                            subject=create_subject_objects(subject_info),
                            sample=create_sample_processing_list(subject_info),
                            data_processing=create_data_processing_list(subject_info)
                            )
    return rep_object


def create_data_processing_list(subject_info):
    """Create a list of DataProcessing objects from subject information."""
    subject_info = fill_missing_required_fields(DataProcessing,  subject_info)

    data_processing_list = []
    data_processing_obj = DataProcessing(data_processing_id=subject_info["data_processing_id"],
                                         # primary_annotation is not required but is not nullable ... feels like a bug in the schema
                                         primary_annotation=subject_info["primary_annotation"] if subject_info["primary_annotation"] else False,
                                         software_versions=subject_info["software_versions"],
                                         paired_reads_assembly=subject_info["paired_reads_assembly"],
                                         quality_thresholds=subject_info["quality_thresholds"],
                                         primer_match_cutoffs=subject_info["primer_match_cutoffs"],
                                         collapsing_method=subject_info["collapsing_method"],
                                         data_processing_protocols=subject_info["data_processing_protocols"],
                                         data_processing_files=[subject_info["data_processing_files"]] if subject_info["data_processing_files"] else None,
                                         germline_database=subject_info["germline_database"],
                                         germline_set_ref=subject_info["germline_set_ref"],
                                         analysis_provenance_id=subject_info["analysis_provenance_id"])
    
    data_processing_list.append(data_processing_obj)
    return data_processing_list
    

def create_sample_processing_list(subject_info):
    """Create a list of SampleProcessing objects from subject information."""
    subject_info = fill_missing_required_fields(Sample,  subject_info)
    subject_info = fill_missing_required_fields(CellProcessing,  subject_info)
    subject_info = fill_missing_required_fields(NucleicAcidProcessing,  subject_info)
    subject_info = fill_missing_required_fields(SequencingRun,  subject_info)

    try:
        library_generation_method = LibraryGenerationMethod(subject_info["library_generation_method"])
    except Exception:
        library_generation_method = LibraryGenerationMethod("other")

    sample_processing_list = []
    sample_processing_obj = SampleProcessing(sample_processing_id=subject_info["sample_processing_id"],
                                             sample_id=subject_info["sample_id"],
                                             sample_type=subject_info["sample_type"],
                                             tissue=Ontology(id=subject_info["tissue_id"], label=subject_info["tissue_label"]),
                                             anatomic_site=subject_info["anatomic_site"],
                                             disease_state_sample=subject_info["disease_state_sample"],
                                             collection_time_point_relative=float(subject_info["collection_time_point_relative"]) if subject_info["collection_time_point_relative"] else 0,
                                             collection_time_point_relative_unit=Ontology(id=subject_info["collection_time_point_relative_unit_id"], label=subject_info["collection_time_point_relative_unit_label"]),
                                             collection_time_point_reference=subject_info["collection_time_point_reference"],
                                             biomaterial_provider=subject_info["biomaterial_provider"],
                                             tissue_processing=subject_info["tissue_processing"],
                                             cell_subset=Ontology(id=subject_info["cell_subset_id"], lable=subject_info["cell_subset_label"]),
                                             cell_phenotype=subject_info["cell_phenotype"],
                                             cell_species=Ontology(id=subject_info["cell_species_id"], lable=subject_info["cell_species_label"]),
                                             single_cell=str_to_bool((subject_info["single_cell"])),
                                             cell_number=subject_info["cell_number"] if subject_info["cell_number"] else 0,
                                             cells_per_reaction=subject_info["cells_per_reaction"] if subject_info["cells_per_reaction"] else 0,
                                             cell_storage=str_to_bool(subject_info["cell_storage"]),
                                             cell_quality=subject_info["cell_quality"],
                                             cell_isolation=subject_info["cell_isolation"],
                                             cell_processing_protocol=subject_info["cell_processing_protocol"],
                                             template_class=TemplateClass(subject_info["template_class"].upper()),
                                             template_quality=subject_info["template_quality"],
                                             template_amount=float(subject_info["template_amount"]) if subject_info["template_amount"] else 0,
                                             template_amount_unit=Ontology(id=subject_info["template_amount_unit_id"], label=subject_info["template_amount_unit_label"]),
                                             library_generation_method=library_generation_method,
                                             library_generation_protocol=subject_info["library_generation_protocol"],
                                             library_generation_kit_version=subject_info["library_generation_kit_version"],
                                             pcr_target=build_pcr_target(subject_info),
                                             complete_sequences=create_complete_sequences_enum(subject_info),
                                             physical_linkage=PhysicalLinkage(subject_info["physical_linkage"]),
                                             sequencing_run_id=subject_info["sequencing_run_id"],
                                             total_reads_passing_qc_filter=subject_info["total_reads_passing_qc_filter"],
                                             sequencing_platform=subject_info["sequencing_platform"],
                                             sequencing_facility=subject_info["sequencing_facility"],
                                             sequencing_run_date=create_date(subject_info, "sequencing_run_date"),
                                             sequencing_kit=subject_info["sequencing_kit"],
                                             sequencing_files=create_sequencing_data_object(subject_info)
                                             )
    
    sample_processing_list.append(sample_processing_obj)

    return sample_processing_list


def build_pcr_target(subject_info):
    targets = []
    for locus in subject_info['pcr_target_locus'].split(','):
        if locus.strip().upper() in ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRG', 'TRD']:
            targets.append(PCRTarget(pcr_target_locus=locus.strip().upper(), 
                                     forward_pcr_primer_target_location=subject_info['forward_pcr_primer_target_location'], 
                                     reverse_pcr_primer_target_location=subject_info['reverse_pcr_primer_target_location']))
    if not targets:
        targets = [PCRTarget(pcr_target_locus=None,
                             forward_pcr_primer_target_location=subject_info['forward_pcr_primer_target_location'], 
                             reverse_pcr_primer_target_location=subject_info['reverse_pcr_primer_target_location'])]
    return targets


def create_complete_sequences_enum(subject_info):
    """Create a CompleteSequences enum from subject information."""
    try:
        return CompleteSequences(subject_info["complete_sequences"])
    except Exception:
        return CompleteSequences("partial")


def create_sequencing_data_object(subject_info):
    """Create a SequencingData object from subject information."""
    subject_info = fill_missing_required_fields(SequencingData,  subject_info)

    sequencing_data_obj = SequencingData(sequencing_data_id=None,
                                         file_type=FileTypeEnum(subject_info["file_type"]) if subject_info["file_type"] else None,
                                         filename=subject_info["filename"],
                                         read_direction=ReadDirectionEnum(subject_info["read_direction"]) if subject_info["read_direction"] else None,
                                         read_length=int(subject_info["read_length"]) if subject_info["read_length"] else None,
                                         paired_filename=subject_info["paired_filename"],
                                         paired_read_direction=PairedReadDirectionEnum(subject_info["paired_read_direction"]) if subject_info["paired_read_direction"] else None,
                                         paired_read_length=int(subject_info["paired_read_length"]) if subject_info["paired_read_length"] else None,
                                         index_filename=None,
                                         index_length=None,
                                         )
    
    return sequencing_data_obj


def str_to_bool(value):
    """Convert a string to a boolean."""
    if isinstance(value, str):
        return value.lower() == 'true'
    return bool(value)


def create_subject_objects(subject_info):
    """Create Subject objects from subject information."""
    subject_info = fill_missing_required_fields(Subject,  subject_info)

    subject_object = Subject(subject_id=subject_info["subject_id"],
                             synthetic=subject_info["synthetic"],
                             species=Ontology(id=subject_info["species_id"], label=subject_info["species_label"]),
                             organism=Ontology(id=subject_info["organism_id"], label=subject_info["organism_label"]), 
                             sex=create_Sex_Enum(subject_info),
                             age_min=subject_info["age_min"],
                             age_max=subject_info["age_max"],
                             age_unit=Ontology(id=subject_info["age_unit_id"], label=subject_info["age_unit_label"]),
                             age_event=subject_info["age_event"],
                             age=subject_info["age"],
                             ancestry_population=subject_info["ancestry_population"],
                             ethnicity=subject_info["ethnicity"],
                             race=subject_info["race"],
                             strain_name=subject_info["strain_name"],
                             linked_subjects=subject_info["linked_subjects"],
                             link_type=subject_info["link_type"],
                             diagnosis=create_diagnosis_list(subject_info),
                             genotype=create_subject_genotype_obj(subject_info))

    return subject_object


def create_Sex_Enum(subject_info):
    """Create a Sex enum from subject information."""
    try:
        return SexEnum(subject_info["sex"])
    except Exception:
        return None


def create_subject_genotype_obj(subject_info):
    """Create a SubjectGenotype object from subject information."""
    subject_genotype_obj = SubjectGenotype(receptor_genotype_set=GenotypeSet(receptor_genotype_set_id=subject_info["receptor_genotype_set_id"] if subject_info["receptor_genotype_set_id"] is not None else "",
                                                                             genotype_class_list=create_genotype_model_list(subject_info)))
    return subject_genotype_obj


def create_genotype_model_list(subject_info):
    """Create a list of GenotypeModel objects from subject information."""
    # genotype_model_list = []
    # genotype_model_obj = GenotypeModel(receptor_genotype_id=)
    return None


def create_diagnosis_list(subject_info):
    """Create a list of Diagnosis objects from subject information."""
    subject_info = fill_missing_required_fields(Diagnosis,  subject_info)

    diagnosis_list = []
    diagnosis_object = Diagnosis(study_group_description=subject_info["study_group_description"],
                                 disease_diagnosis=Ontology(id=subject_info["disease_diagnosis_id"], label=subject_info["disease_diagnosis_label"]),
                                 disease_length=subject_info["disease_length"],
                                 disease_stage=subject_info["disease_stage"],
                                 prior_therapies=subject_info["prior_therapies"],
                                 immunogen=subject_info["immunogen"],
                                 intervention=subject_info["intervention"],
                                 medical_history=subject_info["medical_history"])
    diagnosis_list.append(diagnosis_object)
    return diagnosis_list


def create_study_object(subject_info):
    """Create a Study object from subject information."""
    subject_info = fill_missing_required_fields(Study,  subject_info)

    study_object = Study(study_id=subject_info["study_id"],
                         study_title=subject_info["study_title"],
                         study_type=Ontology(id=subject_info["study_type_id"], label=subject_info["study_type_label"]),
                         study_description=subject_info["study_description"],
                         inclusion_exclusion_criteria=subject_info["inclusion_exclusion_criteria"],
                         grants=subject_info["grants"],
                         study_contact=subject_info["study_contact"],
                         collected_by=subject_info["collected_by"],
                         lab_name=subject_info["lab_name"],
                         lab_address=subject_info["lab_address"],
                         submitted_by=subject_info["submitted_by"],
                         pub_ids=subject_info["pub_ids"],
                         keywords_study=create_keyword_study_list(subject_info),
                         adc_publish_date=create_date(subject_info, "adc_publish_date"),
                         adc_update_date=create_date(subject_info, "adc_update_date"),
                         )
    return study_object


def create_keyword_study_list(subject_info):
    """Create a list of KeywordsStudyEnum from subject information."""
    return [KeywordsStudyEnum('contains_ig')]


def create_date(subject_info, field):
    """Create a datetime object from subject information and field."""
    if not subject_info[field]:
        return None
    else:
        try:
            dateutil.parser.parse(subject_info[field])
        except (ValueError, TypeError):
            return None

        # check for timezone and remove if present

        if '.' in subject_info[field]:
            return subject_info[field].split('.')[0]
        
        return subject_info[field]


def get_default_value(field_type: Any) -> Any:
    """
    Get the default value for a given field type.
    
    Args:
        field_type: The type of the field.
    
    Returns:
        The default value for the field type.
    """
    try:
        if get_origin(field_type) is Union:
            args = get_args(field_type)
            field_type = args[0] if args[1] is type(None) else args[1]

        if field_type == 'int':
            return 0
        elif field_type == 'float':
            return 0.0
        elif field_type == 'str':
            return ""
        elif field_type == 'bool':
            return False
        elif field_type == 'list':
            return []
        elif field_type == 'dict':
            return {}
        elif field_type == 'datetime':
            return datetime.now()
        elif field_type == 'date':
            return datetime.now()
        elif 'Optional' in field_type:
            pass
        else:
            if issubclass(globals()[field_type], Enum):
                # Get the first value of the Enum
                return next(iter(globals()[field_type])).value
    
    except Exception:
        pass

    return None


def fill_missing_required_fields(model_cls: BaseModel, data: dict) -> dict:
    """
    Fill missing required fields in the given data with default values.
    
    Args:
        model_cls: The Pydantic model class.
        data: The data dictionary.
    
    Returns:
        The data dictionary with missing required fields filled.
    """
    filled_data = data.copy()

    for field_name, field_info in model_cls.__fields__.items():
        if field_name in data:
            if is_required(field_info):
                field_type = model_cls.__annotations__[field_name]
                if data[field_name] is None or (field_type != 'str' and data[field_name] == ''):
                    if field_type != 'List[KeywordsStudyEnum]':
                        default_value = get_default_value(field_type)
                    if default_value is not None:
                        filled_data[field_name] = default_value

    return filled_data


def is_required(field_info: FieldInfo) -> bool:
    """
    Check if a field is required.
    
    Args:
        field_info: The field information.
    
    Returns:
        Boolean indicating if the field is required.
    """
    if 'required=True' in str(field_info):
        return True
    
    return False
