import json
from datetime import datetime
from flask import Blueprint, request, jsonify, Response
from schema.models import Enum, date, Ontology, ErrorResponse, SpeciesResponse, Dataset, DatasetsResponse, SubjectDataset, SubjectDatasetResponse, Genotype, Locus, Sample, \
    SampleMetadataResponse, Repertoire, DataProcessing, SampleProcessing, CellProcessing, NucleicAcidProcessing, SequencingRun, LibraryGenerationMethod, TemplateClass, \
    CompleteSequences, PhysicalLinkage, SequencingData, FileTypeEnum, ReadDirectionEnum, PairedReadDirectionEnum, Subject, SexEnum, SubjectGenotype, GenotypeSet, Diagnosis, \
    Study, KeywordsStudyEnum, GenotypeClassListItem, AllSubjectsGenotypeResponse, AllSamplesMetadataResponse
from pydantic.fields import FieldInfo
from pydantic import BaseModel
from typing import Any, Union, get_args, get_origin
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
        json.dumps(encode_obj(obj)),
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
            if 'description' in ds_name and sp not in species_list:
                species_list.append(sp)
                ontology_obj = Ontology(id=ds_data['taxid'], label=ds_data['binomial'])
                ontology_list.append(ontology_obj)

    species_response_obj = SpeciesResponse(species=ontology_list)

    try:
        return species_response_obj.model_dump_json(), 200
    except Exception as e:
        error_response = ErrorResponse(message=str(e))
        return error_response.model_dump_json(), 500


@api_bp.route('/<type>/datasets/<species>', methods=['GET'])
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
def get_subject_datasets(type, species, dataset):
    species = common_lookup(species)

    if not species:
        error_response = ErrorResponse(message="species not found")
        return error_response.model_dump_json(), 500

    """Get subject datasets for a species and dataset based on type."""
    if type == "genomic":
        try:
            subjects_list = genomic.SubjectsAPI(Resource)
            subjects_list = subjects_list.get(species, dataset)
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
            subjects_list = subjects_list.get(species, dataset)
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

            rep_obj = SampleMetadataResponse(repertoire=create_repertoire_obj(sample_info[0]))
            return custom_jsonify(rep_obj.model_dump()), 200

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

            rep_obj = SampleMetadataResponse(repertoire=create_repertoire_obj(sample_info[0]))
            return custom_jsonify(rep_obj.model_dump()), 200

        except Exception as e:
            error_response = ErrorResponse(message=str(e))
            return error_response.model_dump_json(), 500
    else:
        error_response = ErrorResponse(message=str("type not  exists"))
        return error_response.model_dump_json(), 500


@api_bp.route('/<type>/all_samples_metadata/<species>/<dataset>', methods=['GET'])
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
        error_response = ErrorResponse(message=str("type not  exists"))
        return error_response.model_dump_json(), 500

    repertoires = [create_repertoire_obj(s) for s in sample_info[0]]
    rep_obj = AllSamplesMetadataResponse(repertoire_class_list=repertoires)
    return custom_jsonify(rep_obj.model_dump()), 200


def create_repertoire_obj(subject_info):
    """Create a Repertoire object from subject information."""
    subject_info = fill_missing_required_fields(Repertoire,  subject_info)

    rep_object = Repertoire(repertoire_id=subject_info.get("sample_name"),
                            repertoire_name=subject_info.get("repertoire_name", None),
                            repertoire_description=subject_info.get("repertoire_description", None),
                            study=create_study_object(subject_info),
                            subject=create_subject_objects(subject_info),
                            sample=create_sample_processing_list(subject_info),
                            data_processing=create_data_processing_list(subject_info))
    
    return rep_object


def create_data_processing_list(subject_info):
    """Create a list of DataProcessing objects from subject information."""
    subject_info = fill_missing_required_fields(DataProcessing,  subject_info)

    data_processing_list = []
    data_processing_obj = DataProcessing(data_processing_id=subject_info.get("data_processing_id"),
                                         # primary_annotation is not required but is not nullable ... feels like a bug in the schema
                                         primary_annotation=subject_info.get("primary_annotation") if subject_info.get("primary_annotation") else False,
                                         software_versions=subject_info.get("software_versions"),
                                         paired_reads_assembly=subject_info.get("paired_reads_assembly"),
                                         quality_thresholds=subject_info.get("quality_thresholds"),
                                         primer_match_cutoffs=subject_info.get("primer_match_cutoffs"),
                                         collapsing_method=subject_info.get("collapsing_method"),
                                         data_processing_protocols=subject_info.get("data_processing_protocols"),
                                         data_processing_files=[subject_info.get("data_processing_files")] if subject_info.get("data_processing_files") else None,
                                         germline_database=subject_info.get("germline_database"),
                                         germline_set_ref=subject_info.get("germline_set_ref"),
                                         analysis_provenance_id=subject_info.get("analysis_provenance_id"))
    
    data_processing_list.append(data_processing_obj)
    return data_processing_list
    

def create_sample_processing_list(subject_info):
    """Create a list of SampleProcessing objects from subject information."""
    subject_info = fill_missing_required_fields(Sample,  subject_info)
    subject_info = fill_missing_required_fields(CellProcessing,  subject_info)
    subject_info = fill_missing_required_fields(NucleicAcidProcessing,  subject_info)
    subject_info = fill_missing_required_fields(SequencingRun,  subject_info)

    try:
        library_generation_method = LibraryGenerationMethod(subject_info.get("library_generation_method"))
    except Exception:
        library_generation_method = LibraryGenerationMethod("other")

    sample_processing_list = []
    sample_processing_obj = SampleProcessing(sample_processing_id=subject_info.get("sample_processing_id", None),
                                             sample_id=subject_info.get("sample_id"),
                                             sample_type=subject_info.get("sample_type"),
                                             tissue=Ontology(id=subject_info.get("tissue_id", ""), label=subject_info.get("tissue_label", "")),
                                             anatomic_site=subject_info.get("anatomic_site"),
                                             disease_state_sample=subject_info.get("disease_state_sample"),
                                             collection_time_point_relative=int(subject_info.get("collection_time_point_relative")) if subject_info.get("collection_time_point_relative") else 0,
                                             collection_time_point_relative_unit=Ontology(id=subject_info.get("collection_time_point_relative_unit_id", ""), label=subject_info.get("collection_time_point_relative_unit_label", "")),
                                             collection_time_point_reference=subject_info.get("collection_time_point_reference"),
                                             biomaterial_provider=subject_info.get("biomaterial_provider"),
                                             tissue_processing=subject_info.get("tissue_processing"),
                                             cell_subset=Ontology(id=subject_info.get("cell_subset_id", ""), lable=subject_info.get("cell_subset_label", "")),
                                             cell_phenotype=subject_info.get("cell_phenotype"),
                                             cell_species=Ontology(id=subject_info.get("cell_species_id", ""), lable=subject_info.get("cell_species_label", "")),
                                             single_cell=str_to_bool((subject_info.get("single_cell"))),
                                             cell_number=subject_info.get("cell_number") if subject_info.get("cell_number") else 0,
                                             cells_per_reaction=subject_info.get("cells_per_reaction") if subject_info.get("cells_per_reaction") else 0,
                                             cell_storage=str_to_bool(subject_info.get("cell_storage")),
                                             cell_quality=subject_info.get("cell_quality"),
                                             cell_isolation=subject_info.get("cell_isolation"),
                                             cell_processing_protocol=subject_info.get("cell_processing_protocol"),
                                             template_class=TemplateClass(subject_info.get("template_class").upper()),
                                             template_quality=subject_info.get("template_quality"),
                                             template_amount=int(subject_info.get("template_amount")) if subject_info.get("template_amount") else 0,
                                             template_amount_unit=Ontology(id=subject_info.get("template_amount_unit_id", ""), label=subject_info.get("template_amount_unit_label", "")),
                                             library_generation_method=library_generation_method,
                                             library_generation_protocol=subject_info.get("library_generation_protocol"),
                                             library_generation_kit_version=subject_info.get("library_generation_kit_version"),
                                             pcr_target=subject_info.get("pcr_target"),
                                             complete_sequences=create_complete_sequences_enum(subject_info),
                                             physical_linkage=PhysicalLinkage(subject_info.get("physical_linkage")),
                                             sequencing_run_id=subject_info.get("sequencing_run_id"),
                                             total_reads_passing_qc_filter=subject_info.get("total_reads_passing_qc_filter"),
                                             sequencing_platform=subject_info.get("sequencing_platform"),
                                             sequencing_facility=subject_info.get("sequencing_facility"),
                                             sequencing_run_date=subject_info.get("sequencing_run_date"),
                                             sequencing_kit=subject_info.get("sequencing_kit"),
                                             sequencing_files=create_sequencing_data_object(subject_info)
                                             )
    
    sample_processing_list.append(sample_processing_obj)

    return sample_processing_list


def create_complete_sequences_enum(subject_info):
    """Create a CompleteSequences enum from subject information."""
    try:
        return CompleteSequences(subject_info.get("complete_sequences"))
    except Exception:
        return CompleteSequences("partial")


def create_sequencing_data_object(subject_info):
    """Create a SequencingData object from subject information."""
    subject_info = fill_missing_required_fields(SequencingData,  subject_info)

    sequencing_data_obj = SequencingData(sequencing_data_id=subject_info.get("sequencing_data_id") if subject_info.get("sequencing_data_id") is not None else "",
                                         file_type=FileTypeEnum(subject_info.get("file_type")) if subject_info.get("file_type") else None,
                                         filename=subject_info.get("filename"),
                                         read_direction=ReadDirectionEnum(subject_info.get("read_direction")) if subject_info.get("read_direction") else None,
                                         read_length=int(subject_info.get("read_length")) if subject_info.get("read_length") else None,
                                         paired_filename=subject_info.get("paired_filename"),
                                         paired_read_direction=PairedReadDirectionEnum(subject_info.get("paired_read_direction")) if subject_info.get("paired_read_direction") else None,
                                         paired_read_length=int(subject_info.get("paired_read_length")) if subject_info.get("paired_read_length") else None,
                                         index_filename=subject_info.get("index_filename"),
                                         index_length=subject_info.get("index_length"),
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

    subject_object = Subject(subject_id=subject_info.get("subject_id"),
                             synthetic=subject_info.get("synthetic"),
                             species=Ontology(id=subject_info.get("species_id", ""), label=subject_info.get("species_label", "")),
                             organism=Ontology(id=subject_info.get("organism_id", ""), label=subject_info.get("organism_label", "")), 
                             sex=create_Sex_Enum(subject_info),
                             age_min=subject_info.get("age_min"),
                             age_max=subject_info.get("age_max"),
                             age_unit=Ontology(id=subject_info.get("age_unit_id", ""), label=subject_info.get("age_unit_label", "")),
                             age_event=subject_info.get("age_event"),
                             age=subject_info.get("age"),
                             ancestry_population=subject_info.get("ancestry_population"),
                             ethnicity=subject_info.get("ethnicity"),
                             race=subject_info.get("race", None),
                             strain_name=subject_info.get("strain_name", None),
                             linked_subjects=subject_info.get("linked_subjects", None),
                             link_type=subject_info.get("link_type", None),
                             diagnosis=create_diagnosis_list(subject_info),
                             genotype=create_subject_genotype_obj(subject_info))

    return subject_object


def create_Sex_Enum(subject_info):
    """Create a Sex enum from subject information."""
    try:
        return SexEnum(subject_info.get("sex"))
    except Exception:
        return None


def create_subject_genotype_obj(subject_info):
    """Create a SubjectGenotype object from subject information."""
    subject_genotype_obj = SubjectGenotype(receptor_genotype_set=GenotypeSet(receptor_genotype_set_id=subject_info.get("receptor_genotype_set_id") if subject_info.get("receptor_genotype_set_id") is not None else "",
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
    diagnosis_object = Diagnosis(study_group_description=subject_info.get("study_group_description"),
                                 disease_diagnosis=Ontology(id=subject_info.get("disease_diagnosis_id", ""), lable=subject_info.get("disease_diagnosis_label", "")),
                                 disease_length=subject_info.get("disease_length"),
                                 disease_stage=subject_info.get("disease_stage"),
                                 prior_therapies=subject_info.get("prior_therapies"),
                                 immunogen=subject_info.get("immunogen"),
                                 intervention=subject_info.get("intervention"),
                                 medical_history=subject_info.get("medical_history"))
    diagnosis_list.append(diagnosis_object)
    return diagnosis_list


def create_study_object(subject_info):
    """Create a Study object from subject information."""
    subject_info = fill_missing_required_fields(Study,  subject_info)

    study_object = Study(study_id=subject_info.get("study_id"),
                         study_title=subject_info.get("study_title"),
                         study_type=Ontology(id=subject_info.get("study_type_id"), label=subject_info.get("study_type_lable")),
                         study_description=subject_info.get("study_description", None),
                         inclusion_exclusion_criteria=subject_info.get("inclusion_exclusion_criteria", ""),
                         grants=subject_info.get("grants", ""),
                         study_contact=subject_info.get("study_contact", None),
                         collected_by=subject_info.get("collected_by", ""),
                         lab_name=subject_info.get("lab_name", ""),
                         lab_address=subject_info.get("lab_address", ""),
                         submitted_by=subject_info.get("submitted_by", ""),
                         pub_ids=subject_info.get("pub_ids", ""),
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
    if subject_info.get(field) == "":
        return None
    
    else:
        return datetime.now()


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
